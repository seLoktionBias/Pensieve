#!/usr/bin/env python3
import argparse, csv, difflib, subprocess, sys
from collections import Counter
from pathlib import Path
from Bio import Phylo, SeqIO
from Bio.Align import PairwiseAligner

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


def build_raw_to_alignment_map(seq, raw_seq):
    """Build a full 1-based-raw-position -> 1-based-alignment-column map for
    one species' own row, once, for repeated lookup (see map_raw_nt_to_aln).

    MACSE never deletes or reorders a real input character, but it can
    insert synthetic characters that are not present in the raw input at
    all, in more than one place in the same row -- not only its own
    documented '!' partial-codon marker, but (confirmed on real data:
    GUCA1C/Desmodus_rotundus) occasionally a literal ambiguity code like
    'N' standing in for a position it could not resolve, with no '!'
    involved, and sometimes more than once in the same row. A naive
    single-pass greedy character-by-character match (skip on mismatch,
    otherwise advance) breaks down as soon as there is more than one such
    insertion: a mismatch at the first one is handled correctly, but if the
    very next real raw character coincidentally does not equal whatever
    comes next in seq (extremely likely once there is a second, unrelated
    insertion later in the row), the greedy scan has no way to recover --
    it can stall indefinitely with no further matches, or silently latch
    onto the wrong character. This showed up on real data: with a second,
    independent MACSE insertion further into the same gene, greedy matching
    for Desmodus_rotundus stalled partway through GUCA1C and never found a
    STOP coordinate mapping in the back half of the gene at all.

    A proper sequence alignment (difflib's longest-matching-block algorithm)
    between the species' own raw sequence and the non-gap characters of its
    aligned row does not have this failure mode: every 'equal' opcode block
    is, by construction, a genuine, correctly-ordered run of real raw
    characters, however many separate synthetic insertions sit between them.
    """
    non_gap_positions = [i for i, c in enumerate(seq) if c != "-"]
    non_gap_chars = "".join(seq[i] for i in non_gap_positions)
    sm = difflib.SequenceMatcher(None, raw_seq.upper(), non_gap_chars.upper(), autojunk=False)
    mapping = {}
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag != "equal":
            continue
        for k in range(i2 - i1):
            mapping[i1 + k + 1] = non_gap_positions[j1 + k] + 1
    return mapping


def map_raw_nt_to_aln(seq, raw_nt, raw_seq):
    """Map a single 1-based raw position to its 1-based alignment column.
    Convenience wrapper around build_raw_to_alignment_map() for one-off
    lookups; prefer building the map once and reusing it for many lookups
    against the same (seq, raw_seq) pair.
    """
    return build_raw_to_alignment_map(seq, raw_seq).get(raw_nt)


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
    complete_fasta = step00 / f"00_{gene}.complete_seqs.fa"
    incomplete_fasta = step00 / f"00_{gene}.incomplete_seqs.fa"
    nt_out = outdir / f"01_{gene}.macse_NT.fasta"
    aa_out = outdir / f"01_{gene}.macse_AA.fasta"

    # MACSE is told which sequences it may trust, instead of being second-guessed
    # afterwards.
    #
    # Complete-ORF sequences (step 00's 00_<gene>.complete_seqs.fa) go to -seq,
    # MACSE's "reliable" slot; every incomplete/pseudogenized sequence goes to
    # -seq_lr, which MACSE documents as the slot for pseudogenes. A frameshift
    # inside a reliable sequence then costs FRAMESHIFT_COST_RELIABLE, and one at
    # its first or last codon costs FRAMESHIFT_COST_RELIABLE_TERMINAL: both are
    # set high enough that MACSE will always prefer any other explanation, so an
    # intact reading frame cannot be broken open to make a genuinely broken
    # lineage align more cheaply. Frameshifts in the less-reliable set keep
    # MACSE's own default cost, which is the entire point -- that is where the
    # real lesions are.
    #
    # This replaces the conserved 5' start-block heuristic (removed): that rule
    # detected the same corruption after the fact, but only in the first
    # 30 codons, only for genes with at least five complete ORFs, and only at
    # the 5' end. Declaring reliability up front prevents it everywhere in the
    # gene instead of repairing it in one window.
    FRAMESHIFT_COST_RELIABLE = "1000"
    FRAMESHIFT_COST_RELIABLE_TERMINAL = "1000"

    def run_macse(log_name):
        tool_cwd = outdir / "_external_tool_work" / "macse"
        n_complete = sum(1 for _ in read_fasta(complete_fasta)) if complete_fasta.exists() else 0
        n_incomplete = sum(1 for _ in read_fasta(incomplete_fasta)) if incomplete_fasta.exists() else 0
        cmd = [args.macse_cmd, "-prog", "alignSequences"]
        if n_complete:
            cmd += ["-seq", str(complete_fasta.resolve())]
            if n_incomplete:
                cmd += ["-seq_lr", str(incomplete_fasta.resolve())]
        else:
            # -seq is mandatory and -seq_lr alone is not a valid MACSE input, so
            # a gene with no complete ORF at all falls back to the undivided
            # input with MACSE's default costs. Raising -fs there would forbid
            # frameshifts in sequences that are known to be frameshifted.
            print(f"[WARN] {gene}: no complete-ORF sequence; running MACSE on the undivided input "
                  f"with default frameshift costs.", file=sys.stderr)
            cmd += ["-seq", str(Path(fasta).resolve())]
        if n_complete:
            cmd += ["-fs", FRAMESHIFT_COST_RELIABLE,
                    "-fs_term", FRAMESHIFT_COST_RELIABLE_TERMINAL]
        cmd += ["-out_NT", str(nt_out.resolve()), "-out_AA", str(aa_out.resolve())]
        print(f"Running MACSE ({n_complete} reliable / {n_incomplete} less-reliable sequence(s)); "
              f"log: {outdir / log_name}", file=sys.stderr)
        run_external(cmd, tool_cwd, outdir / log_name)

    if not nt_out.exists() or not aa_out.exists():
        run_macse(f"01_{gene}.macse.log")

    raw = read_fasta(fasta)
    status = read_tsv(step00 / f"00_{gene}.orf_status.tsv")
    status_by = {r["species"]: r for r in status if "species" in r}

    nt = read_fasta(nt_out)
    aa = read_fasta(aa_out)

    # MACSE may preserve, alter, or alphabetize record order depending on version.
    # Rewrite both MACSE outputs in the tree order that matches how
    # 05_plot_events.R actually renders the tree top-to-bottom (see
    # terminals_top_to_bottom() in 00_prune_and_check_orf.py for why plain
    # left-to-right preorder reads bottom-to-top next to the plot instead).
    tree = Phylo.read(str(step00 / f"00_{gene}.common_species.tree"), "newick")
    def terminals_top_to_bottom(clade):
        if not clade.clades:
            return [clade]
        tips = []
        for child in reversed(clade.clades):
            tips.extend(terminals_top_to_bottom(child))
        return tips
    phylo_order = [tip.name for tip in terminals_top_to_bottom(tree.root)]
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
    raw_to_aln_by_species = {}
    for r in failures:
        if r.get("failure_type") != "premature_in_frame_stop":
            continue
        sp = r["species"]
        seq = nt.get(sp)
        if not seq:
            continue
        if sp not in raw_to_aln_by_species:
            raw_to_aln_by_species[sp] = build_raw_to_alignment_map(seq, raw[sp])
        raw_to_aln = raw_to_aln_by_species[sp]
        nt_start = int(float(r["nt_start"]))
        nt_end = int(float(r["nt_end"]))
        cols = [raw_to_aln.get(x) for x in range(nt_start, nt_end + 1)]
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
        f"Finished step01 for {gene}; species={len(nt)}; alignment_length={len(next(iter(nt.values())))}; "
        f"frameshift_markers={len(fs)}"
    )


if __name__ == "__main__":
    main()
