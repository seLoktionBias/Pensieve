#!/usr/bin/env python3
import argparse, csv, subprocess, sys
from pathlib import Path
from Bio import Phylo, SeqIO

STOP = {"TAA", "TAG", "TGA"}


def read_tsv(path):
    if not Path(path).exists():
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(rows, path, header=None):
    if header is None:
        header = []
        for r in rows:
            for k in r:
                if k not in header:
                    header.append(k)
    with open(path, "w", newline="") as out:
        w = csv.DictWriter(out, delimiter="\t", fieldnames=header, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_fasta(path):
    d = {}
    for r in SeqIO.parse(path, "fasta"):
        d[r.id.split()[0]] = str(r.seq).upper()
    return d


def write_fasta(records, path, order, width=80):
    with open(path, "w") as out:
        for name in order:
            if name not in records:
                continue
            out.write(f">{name}\n")
            seq = records[name]
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def contiguous_runs(seq, chars):
    runs = []
    i = 0
    n = len(seq)
    chars = set(chars)
    while i < n:
        if seq[i] in chars:
            j = i
            while j < n and seq[j] in chars:
                j += 1
            runs.append((i + 1, j, seq[i:j]))
            i = j
        else:
            i += 1
    return runs


def species_raw_pos_at(seq, aln_start, aln_end):
    """Map an alignment interval to flanking/raw coordinates for that same species.

    This is purely within-species provenance. It never compares against another
    sequence and therefore does not create a reference taxon or reference-relative
    insertion/deletion label.
    """
    left = sum(1 for c in seq[:aln_start - 1] if c not in "-!")
    within = sum(1 for c in seq[aln_start - 1:aln_end] if c not in "-!")
    if within == 0:
        return left, left
    return left + 1, left + within


def map_raw_nt_to_aln(seq, raw_nt):
    count = 0
    for i, c in enumerate(seq, start=1):
        if c in "-!":
            continue
        count += 1
        if count == raw_nt:
            return i
    return None


def run_external(cmd, cwd, log_file):
    cwd = Path(cwd)
    cwd.mkdir(parents=True, exist_ok=True)
    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with open(log_file, "w") as log:
        print("Running:", " ".join(str(x) for x in cmd), file=log)
        print("Working directory:", str(cwd), file=log)
        log.flush()
        try:
            subprocess.check_call(cmd, cwd=str(cwd), stdout=log, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] External command failed with exit code {e.returncode}. See log: {log_file}", file=sys.stderr)
            raise


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--step00-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--macse-cmd", default="macse")
    args = ap.parse_args()

    gene = args.gene
    step00 = Path(args.step00_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    macse_input = step00 / f"00_{gene}.common_species.gapless_for_macse.fasta"
    fasta = macse_input if macse_input.exists() else step00 / f"00_{gene}.common_species.fasta"
    nt_out = outdir / f"01_{gene}.macse_NT.fasta"
    aa_out = outdir / f"01_{gene}.macse_AA.fasta"

    if not nt_out.exists() or not aa_out.exists():
        tool_cwd = outdir / "_external_tool_work" / "macse"
        cmd = [
            args.macse_cmd, "-prog", "alignSequences", "-seq", str(Path(fasta).resolve()),
            "-out_NT", str(nt_out.resolve()), "-out_AA", str(aa_out.resolve())
        ]
        print("Running MACSE; log:", outdir / f"01_{gene}.macse.log", file=sys.stderr)
        run_external(cmd, tool_cwd, outdir / f"01_{gene}.macse.log")

    raw = read_fasta(fasta)
    nt = read_fasta(nt_out)
    aa = read_fasta(aa_out)

    # MACSE may preserve, alter, or alphabetize record order depending on version.
    # Rewrite both MACSE outputs in the left-to-right order of the authoritative
    # pruned species tree so every inspectable MSA follows phylogeny.
    tree = Phylo.read(str(step00 / f"00_{gene}.common_species.tree"), "newick")
    phylo_order = [tip.name for tip in tree.get_terminals()]
    if set(phylo_order) != set(nt):
        raise SystemExit(
            f"MACSE NT species differ from pruned tree: missing={sorted(set(phylo_order)-set(nt))[:10]}, "
            f"extra={sorted(set(nt)-set(phylo_order))[:10]}"
        )
    if set(aa) != set(nt):
        raise SystemExit("MACSE AA and NT outputs contain different species")
    write_fasta(nt, nt_out, phylo_order)
    write_fasta(aa, aa_out, phylo_order)
    nt = {name: nt[name] for name in phylo_order}
    aa = {name: aa[name] for name in phylo_order}

    if not nt:
        raise SystemExit(f"MACSE NT alignment is empty: {nt_out}")
    lengths = {len(seq) for seq in nt.values()}
    if len(lengths) != 1:
        raise SystemExit(f"MACSE NT output is not a rectangular alignment: lengths={sorted(lengths)}")

    status = read_tsv(step00 / f"00_{gene}.orf_status.tsv")
    status_by = {r["species"]: r for r in status if "species" in r}
    rows = []
    for sp, seq in nt.items():
        rawseq = raw.get(sp, "")
        rows.append({
            "gene": gene,
            "species": sp,
            "raw_length": len(rawseq),
            "macse_nt_length": len(seq),
            "macse_gap_count": seq.count("-"),
            "macse_frameshift_marker_count": seq.count("!"),
            "macse_aa_stop_count": aa.get(sp, "").count("*"),
            "raw_premature_stop_count": status_by.get(sp, {}).get("premature_stop_count", "NA"),
            "raw_complete_orf": status_by.get(sp, {}).get("complete_orf", "NA"),
            "coordinate_system": "macse_alignment_1_based_inclusive",
        })
    write_tsv(rows, outdir / f"01_{gene}.macse_sequence_summary.tsv")

    # MACSE '!' is a partial-codon / frameshift placeholder. We record its
    # location in MACSE alignment coordinates without assigning insertion or
    # deletion polarity and without comparing the sequence to any reference.
    fs = []
    for sp, seq in nt.items():
        for s, e, seg in contiguous_runs(seq, "!"):
            ss, se = species_raw_pos_at(seq, s, e)
            fs.append({
                "gene": gene,
                "species": sp,
                "event_type": "macse_partial_codon_marker",
                "alignment_start": s,
                "alignment_end": e,
                "marker_length": e - s + 1,
                "species_raw_left_nt": ss,
                "species_raw_right_nt": se,
                "coordinate_system": "macse_alignment_1_based_inclusive",
                "classification": "frameshift_placeholder_no_polarity",
            })
    write_tsv(
        fs,
        outdir / f"01_{gene}.macse_frameshift_markers.tsv",
        ["gene", "species", "event_type", "alignment_start", "alignment_end", "marker_length",
         "species_raw_left_nt", "species_raw_right_nt", "coordinate_system", "classification"],
    )

    # Raw premature STOPs are mapped onto MACSE alignment coordinates before
    # PAML masking. They remain allele-specific downstream. No reference taxon
    # or reference-coordinate projection is involved.
    failures = read_tsv(step00 / f"00_{gene}.orf_failures.tsv")
    pm = []
    for r in failures:
        if r.get("failure_type") != "premature_in_frame_stop":
            continue
        sp = r["species"]
        seq = nt.get(sp)
        if not seq:
            continue
        nt_start = int(float(r["nt_start"]))
        nt_end = int(float(r["nt_end"]))
        cols = [map_raw_nt_to_aln(seq, x) for x in range(nt_start, nt_end + 1)]
        if any(c is None for c in cols):
            raise SystemExit(
                f"Could not map raw premature STOP coordinates for {sp}: {nt_start}-{nt_end} into {nt_out}"
            )
        upstream = sorted(
            [x for x in fs if x["species"] == sp and int(x["alignment_end"]) < min(cols)],
            key=lambda x: (int(x["alignment_end"]), int(x["alignment_start"])),
        )
        nearest = upstream[-1] if upstream else None
        # MACSE ! runs encode the number of artificial nucleotides needed to
        # restore codon phase in the alignment.  Their cumulative length modulo
        # three therefore tells us whether the reading frame is shifted at this
        # STOP.  Crucially, "any upstream !" is NOT enough: a later compensating
        # frameshift can restore phase and a downstream nonsense mutation must
        # then remain an independent STOP candidate.
        frame_correction_mod3 = sum(int(x["marker_length"]) for x in upstream) % 3
        frame_shifted_at_stop = frame_correction_mod3 != 0
        if frame_shifted_at_stop:
            interpretation = "raw_stop_likely_frameshift_consequence"
        elif upstream:
            interpretation = "raw_stop_independent_candidate_after_compensated_frame_restoration"
        else:
            interpretation = "raw_stop_independent_candidate"
        pm.append({
            "gene": gene,
            "species": sp,
            "codon_position": r.get("codon_position"),
            "raw_nt_start": nt_start,
            "raw_nt_end": nt_end,
            "raw_stop_codon": r.get("codon"),
            "alignment_columns": ",".join(str(c) for c in cols),
            "coordinate_system": "macse_alignment_1_based_inclusive",
            "upstream_macse_marker_count": len(upstream),
            "upstream_macse_frame_correction_mod3": frame_correction_mod3,
            "frame_shifted_at_stop": frame_shifted_at_stop,
            # Compatibility field retained for old workdirs/readers. Its v3.30
            # meaning is now "frame is shifted at STOP", not merely "a marker
            # exists somewhere upstream".
            "masked_by_upstream_macse_frameshift_marker": frame_shifted_at_stop,
            "upstream_marker_alignment_start": nearest["alignment_start"] if nearest else "NA",
            "upstream_marker_alignment_end": nearest["alignment_end"] if nearest else "NA",
            "interpretation": interpretation,
        })
    write_tsv(
        pm,
        outdir / f"01_{gene}.macse_premature_stop_masking.tsv",
        ["gene", "species", "codon_position", "raw_nt_start", "raw_nt_end", "raw_stop_codon",
         "alignment_columns", "coordinate_system", "upstream_macse_marker_count",
         "upstream_macse_frame_correction_mod3", "frame_shifted_at_stop",
         "masked_by_upstream_macse_frameshift_marker", "upstream_marker_alignment_start",
         "upstream_marker_alignment_end", "interpretation"],
    )

    print(
        f"Finished step01 for {gene}; species={len(nt)}; alignment_length={next(iter(lengths))}; "
        f"frameshift_markers={len(fs)}; coordinate_system=MACSE alignment"
    )


if __name__ == "__main__":
    main()
