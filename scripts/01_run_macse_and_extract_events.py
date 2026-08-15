#!/usr/bin/env python3
import argparse, csv, subprocess, sys
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


def map_raw_nt_to_aln(seq, raw_nt):
    count = 0
    for i, c in enumerate(seq, start=1):
        if c in "-!":
            continue
        count += 1
        if count == raw_nt:
            return i
    return None


STANDARD_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_codon_seq(seq):
    """Translate a gapless codon sequence to amino acids. Any codon
    containing a non-ACGT character (in particular the N used to mask an
    untrusted MACSE modification) translates to the ambiguous residue 'X'
    rather than guessing a specific amino acid from partial information.
    """
    aa = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        aa.append(STANDARD_CODON_TABLE.get(seq[i:i + 3], "X"))
    return "".join(aa)


def is_complete_orf(status_by, name):
    return str(status_by.get(name, {}).get("complete_orf", "")).strip().lower() in ("true", "1")


def detect_conserved_block(raw, status_by, max_probe_codons=30, min_species=5, identity_threshold=0.3):
    """How many codons, counting from position 1, do the complete-ORF
    sequences agree on closely enough to call "the conserved start block" --
    the region MACSE most often corrupts by misplacing a frame-correction
    marker meant for one pseudogenized lineage's real, more downstream
    lesion (see resolve_conserved_block). Detection uses ONLY complete-ORF
    raw sequences, which by construction all start at their own real codon 1
    with no gaps needed between them this early in the gene, so no
    alignment step is required to compare them position by position.
    Returns 0 (block correction disabled) if there are too few complete
    ORFs to establish a reliable consensus.

    identity_threshold is deliberately low (real PDE6C data: modal amino
    acid identity across complete ORFs bounces around at ordinary,
    non-indel polymorphic sites -- e.g. a genuine Ile/Val substitution
    position at only 60% identity, sandwiched between 95-100% identity
    neighbours); the identity check exists only to stop extending the block
    into a region so divergent it suggests a genuine indel breaking
    positional homology between species (never observed above ~0.29 in real
    gene data probed so far), not to demand near-total conservation at
    every single codon.

    "complete ORF" (starts with ATG, length a multiple of 3, no internal
    in-frame STOP -- see 00_prune_and_check_orf.py) is an ORF-SCREEN result,
    not a structural-alignment proof: a compensating pair of frameshifts
    could in principle still satisfy it. It is used here only as a cheap,
    conservative filter for "safe enough to compare position-by-position
    without an alignment step first," not as a guarantee that these
    sequences are internally collinear at every position.
    """
    complete = [s for n, s in raw.items() if is_complete_orf(status_by, n)]
    if len(complete) < min_species:
        return 0
    max_codons = min(max_probe_codons, min(len(s) // 3 for s in complete))
    block_codons = 0
    for i in range(max_codons):
        aas = [translate_codon_seq(s[i * 3:i * 3 + 3]) for s in complete]
        modal_aa, modal_count = Counter(aas).most_common(1)[0]
        if modal_aa == "X" or modal_count / len(aas) < identity_threshold:
            break
        block_codons += 1
    return block_codons


def block_reference(raw, status_by, block_codons):
    """Majority-vote codon at each of the first block_codons positions,
    built only from complete-ORF sequences (see detect_conserved_block).

    This is a symmetric, computed STATISTIC over every complete-ORF
    species treated equally -- not a reference SPECIES: no single taxon's
    own sequence is privileged, no coordinate system is projected from it,
    and it plays no role in event polarity/direction (Sankoff parsimony on
    the phylogeny decides that, downstream and independently). It exists
    purely as an internal QC baseline for block_quality_scores() to compare
    MACSE's proposed edit against, in resolve_conserved_block() below.
    """
    complete = [s for n, s in raw.items() if is_complete_orf(status_by, n)]
    ref_codons = []
    for i in range(block_codons):
        codons = [s[i * 3:i * 3 + 3] for s in complete]
        ref_codons.append(Counter(codons).most_common(1)[0][0])
    return "".join(ref_codons)


_BLOCK_ALIGNER = PairwiseAligner()
_BLOCK_ALIGNER.mode = "global"
_BLOCK_ALIGNER.match_score = 1
_BLOCK_ALIGNER.mismatch_score = 0
_BLOCK_ALIGNER.open_gap_score = -1
_BLOCK_ALIGNER.extend_gap_score = -0.5


def block_quality_scores(candidate_nt, ref_nt):
    """Normalised global-alignment identity of a candidate block sequence
    against the conserved-block reference, at both the nucleotide and
    amino-acid level (both roughly 0..1; higher is better)."""
    nt_score = _BLOCK_ALIGNER.score(candidate_nt, ref_nt) / max(len(ref_nt), 1)
    aa_score = (_BLOCK_ALIGNER.score(translate_codon_seq(candidate_nt), translate_codon_seq(ref_nt))
                / max(len(ref_nt) // 3, 1))
    return nt_score, aa_score


def resolve_conserved_block(raw, nt_macse0, status_by, block_codons):
    """Decide, for every species, what its own conserved-block content
    should be: its own untouched raw block prefix by default, or MACSE's
    own proposed edit there, ONLY if that edit demonstrably improves
    alignment quality against the block reference.

    Real bug, found by inspecting real PDE6C data directly: MACSE can place
    its frame-restoration marker anywhere a global DP alignment finds it
    cheapest, not necessarily near a lineage's actual lesion. Several
    vampire-bat-clade species (e.g. Desmodus_rotundus, Diaemus_youngii)
    whose own raw sequence starts with a completely normal ATG -- their real
    pseudogenizing lesion is elsewhere in the gene -- nonetheless got an
    extra frame-correction codon inserted right at the conserved 5' start by
    MACSE, corrupting a region that never needed correcting; blindly
    trusting MACSE's own row for every non-complete-ORF tip propagated that
    corruption downstream, spreading it to every OTHER species' start codon
    too via the shared alignment columns (the original PDE6C "A--TG" bug).
    Simply re-running MACSE on a "corrected" version of the full input does
    NOT fix this by itself: MACSE sees the exact same conserved region again
    and reliably makes the exact same placement decision again. The block
    region must never be shown to MACSE at all for it to stop being
    corrupted (see main(), which strips the decided block content off the
    front of every sequence before MACSE ever runs on the remainder).

    For each species with a MACSE modification inside the block: build two
    candidates -- the species' own untouched raw block prefix, and the
    portion of MACSE's own alignment row spanning that same block (degapped,
    '!' -> N) -- and score both against the block reference consensus (see
    block_reference), on nucleotide AND amino-acid identity. MACSE's
    modification is kept only if it does not score worse on either axis and
    scores strictly better on at least one; otherwise the species' own raw
    block content is used.

    Returns (block_content, log_rows). block_content maps every species to
    its own decided block string (NOT necessarily block_nt long -- an
    accepted MACSE edit may differ in length from a genuine indel; see
    realign_block_content()).
    """
    block_nt = block_codons * 3
    ref = block_reference(raw, status_by, block_codons)
    block_content = {}
    log_rows = []
    for name, raw_seq in raw.items():
        block_content[name] = raw_seq[:block_nt] if len(raw_seq) >= block_nt else raw_seq
        macse_seq = nt_macse0.get(name)
        if macse_seq is None or len(raw_seq) < block_nt:
            continue
        aln_col_end = map_raw_nt_to_aln(macse_seq, block_nt)
        if aln_col_end is None:
            continue
        block_region = macse_seq[:aln_col_end]
        if "!" not in block_region and "-" not in block_region:
            continue  # MACSE made no modification in the block for this species
        candidate_raw = raw_seq[:block_nt]
        candidate_macse = block_region.replace("!", "N").replace("-", "")
        raw_nt_score, raw_aa_score = block_quality_scores(candidate_raw, ref)
        macse_nt_score, macse_aa_score = block_quality_scores(candidate_macse, ref)
        accept = (macse_nt_score >= raw_nt_score and macse_aa_score >= raw_aa_score
                  and (macse_nt_score > raw_nt_score or macse_aa_score > raw_aa_score))
        if accept:
            block_content[name] = candidate_macse
        log_rows.append({
            "species": name, "block_codons": block_codons, "macse_modification_accepted": accept,
            "raw_block_nt_identity": round(raw_nt_score, 3), "raw_block_aa_identity": round(raw_aa_score, 3),
            "macse_block_nt_identity": round(macse_nt_score, 3), "macse_block_aa_identity": round(macse_aa_score, 3),
        })
    return block_content, log_rows


def realign_block_content(block_content, macse_cmd, tool_cwd, log_file):
    """block_content values are usually all exactly block_nt long (every
    species either kept its own raw block prefix or had no MACSE
    modification there at all) -- already a valid, gapless mini-alignment
    by construction, so this returns immediately with no extra work in that
    (overwhelmingly common) case.

    Only when an ACCEPTED MACSE edit genuinely changed a species' block
    length (a real indel inside the conserved block, or its own count of
    MACSE '!' frame-restoration characters folded into that length) does
    this need to actually realign. Pensieve v4.2 uses MACSE itself for this
    -- the single alignment engine everywhere in the pipeline, never a
    second tool -- instead of the earlier approach of translating to
    protein and aligning with a separate amino-acid aligner. MACSE aligns
    raw nucleotide block candidates directly and is frame-aware by
    construction, so it does not have the earlier bug (found running real
    CNGA3 data) of producing a block width not divisible by 3 when
    candidates differ in length. Each candidate is first trimmed to its own largest multiple of
    3 (an accepted MACSE block edit is already codon-respecting, so this
    can only discard 1-2 trailing nucleotides of an already-provisional
    block candidate, never real, already-accepted sequence content).
    """
    lengths = {len(seq) for seq in block_content.values()}
    if len(lengths) == 1:
        return dict(block_content)
    trimmed = {name: seq[:len(seq) - len(seq) % 3] for name, seq in block_content.items()}

    tool_cwd = Path(tool_cwd); tool_cwd.mkdir(parents=True, exist_ok=True)
    macse_in = tool_cwd / "block_input.fasta"
    nt_out = tool_cwd / "block_aligned_NT.fasta"
    aa_out = tool_cwd / "block_aligned_AA.fasta"
    write_fasta(trimmed, macse_in, list(trimmed.keys()))
    cmd = [macse_cmd, "-prog", "alignSequences", "-seq", str(macse_in.resolve()),
           "-out_NT", str(nt_out.resolve()), "-out_AA", str(aa_out.resolve())]
    run_external(cmd, tool_cwd, log_file)
    nt_aligned = read_fasta(nt_out)
    if set(nt_aligned) != set(block_content):
        raise SystemExit("MACSE block realignment species mismatch")
    aligned_lengths = {len(seq) for seq in nt_aligned.values()}
    if len(aligned_lengths) != 1 or next(iter(aligned_lengths)) % 3 != 0:
        raise SystemExit(f"MACSE block realignment is not a codon-multiple rectangular alignment: "
                          f"lengths={sorted(aligned_lengths)}")
    return nt_aligned


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

    def run_macse(seq_fasta, log_name):
        tool_cwd = outdir / "_external_tool_work" / "macse"
        cmd = [
            args.macse_cmd, "-prog", "alignSequences", "-seq", str(Path(seq_fasta).resolve()),
            "-out_NT", str(nt_out.resolve()), "-out_AA", str(aa_out.resolve())
        ]
        print("Running MACSE; log:", outdir / log_name, file=sys.stderr)
        run_external(cmd, tool_cwd, outdir / log_name)

    if not nt_out.exists() or not aa_out.exists():
        run_macse(fasta, f"01_{gene}.macse.log")

    raw = read_fasta(fasta)
    status = read_tsv(step00 / f"00_{gene}.orf_status.tsv")
    status_by = {r["species"]: r for r in status if "species" in r}

    # Conserved 5' block quality-gate: MACSE can place its frame-correction
    # marker anywhere a global DP alignment finds cheapest, not necessarily
    # near a lineage's real lesion, and has been observed (real PDE6C data)
    # corrupting the highly conserved start block for species whose actual
    # pseudogenizing lesion lies elsewhere in the gene. See
    # resolve_conserved_block() for the full rationale on WHY the block is
    # spliced out and MACSE re-run on the remainder only, instead of simply
    # re-running MACSE on a "corrected" full sequence (which does not work:
    # MACSE sees the identical conserved region again and reliably makes the
    # identical, wrong placement decision again).
    nt0 = read_fasta(nt_out)
    block_codons = detect_conserved_block(raw, status_by)
    block_content, block_log = resolve_conserved_block(raw, nt0, status_by, block_codons)
    write_tsv(
        block_log, outdir / f"01_{gene}.conserved_block_correction.tsv",
        ["species", "block_codons", "macse_modification_accepted",
         "raw_block_nt_identity", "raw_block_aa_identity",
         "macse_block_nt_identity", "macse_block_aa_identity"],
    )
    if block_codons >= 3:
        block_content = realign_block_content(
            block_content, args.macse_cmd,
            outdir / "_external_tool_work" / "block_realign", outdir / f"01_{gene}.block_realign.log",
        )
        block_nt = block_codons * 3
        remainder_raw = {name: raw[name][block_nt:] for name in raw}
        remainder_fasta = outdir / f"01_{gene}.block_remainder_input.fasta"
        write_fasta(remainder_raw, remainder_fasta, list(raw.keys()))
        n_flagged = sum(1 for r in block_log if r["macse_modification_accepted"] in (True, "True"))
        print(f"Conserved block check for {gene}: first {block_codons} codons are never shown to "
              f"MACSE (each species gets its own quality-gated block content directly; "
              f"{n_flagged} MACSE block edits were accepted, {len(block_log) - n_flagged} rejected "
              f"in favour of the species' own raw sequence); MACSE runs only on the remainder.",
              file=sys.stderr)
        run_macse(remainder_fasta, f"01_{gene}.macse_remainder.log")
        remainder_nt = read_fasta(nt_out)
        remainder_aa = read_fasta(aa_out)
        nt = {name: block_content[name] + remainder_nt[name] for name in remainder_nt}
        aa = {name: translate_codon_seq(block_content[name]).replace("X", "-") + remainder_aa[name]
              for name in remainder_aa}
    else:
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
        f"Finished step01 for {gene}; species={len(nt)}; alignment_length={len(next(iter(nt.values())))}; "
        f"frameshift_markers={len(fs)}"
    )


if __name__ == "__main__":
    main()
