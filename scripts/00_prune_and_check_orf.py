#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path
from Bio import Phylo, SeqIO

STOP_CODONS = {"TAA", "TAG", "TGA"}


def write_tsv(rows, path, header):
    with open(path, "w", newline="") as out:
        w = csv.DictWriter(out, delimiter="\t", fieldnames=header, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def clean_seq(seq):
    """Keep biological sequence characters and user-defined alignment gaps.

    N is retained because users often use it for assembly uncertainty.  Gaps
    are retained in the stored FASTA so --alignment defined can preserve a
    manually curated alignment, but ORF checks below remove gaps first.
    """
    seq = str(seq).upper().replace("U", "T")
    return "".join(c for c in seq if c in "ACGTN-")


def gapless_seq(seq):
    return clean_seq(seq).replace("-", "")


def strip_terminal_stop(seq, preserve_width=False):
    """Remove a single trailing TAA/TAG/TGA -- whatever this sequence's own
    real (gap-stripped) last 3 characters are, if they spell a stop codon --
    before anything else in the pipeline ever sees it.

    preserve_width (used for --alignment defined): the input is one row of a
    user-curated codon alignment whose column count is authoritative and must
    stay identical across every row. Deleting the 3 terminal characters there
    would shorten only the rows that happen to carry a terminal stop -- and,
    because a defined row's real terminal codon is frequently followed by
    trailing alignment gaps, would delete characters from the MIDDLE of the row
    rather than its end -- desynchronising the alignment and making step 02's
    equal-length check fail (confirmed on real PDE6C/CNGB3/GUCY2F data, where
    only some species carry a terminal stop). With preserve_width the terminal
    stop's three real characters are replaced by gaps in place instead, so the
    stop is genuinely removed from every downstream (gap-stripped) view while
    the alignment keeps exactly the same columns Pensieve promised never to
    add or remove for a defined alignment.

    v4.0: Pensieve no longer reconstructs a terminal stop codon at ancestral
    nodes at all (see CHANGELOG) -- PAML's codon model structurally cannot
    represent a stop as a state, and every attempt at parsimony-reconstructing
    it from the tip alignment (three different designs across v3.35-v3.37)
    kept surfacing new real failure modes as the alignment coordinate system
    scattered a tip's own terminal codon in ways each fix only partially
    covered. Removing the terminal stop from every input sequence up front
    sidesteps the whole problem instead of continuing to chase it: there is
    nothing left for any downstream step to reconstruct, mis-scatter, or
    disagree about at the very end of the CDS.

    Deliberately does NOT require len(real) % 3 == 0: a pseudogenized
    sequence carrying an internal frameshift has no single consistent
    reading frame end-to-end, but its literal last 3 real characters can
    still spell a stop codon that must be stripped like any other -- gating
    on mod3 left frameshifted sequences (e.g. Desmodus_rotundus, real PDE6C
    data) with an untouched trailing TAA.

    Removes exactly the sequence's own last 3 real (non-gap) characters,
    wherever they fall in the string -- not just the literal last 3
    characters -- so this is safe for both a gapless raw CDS and a
    user-curated --alignment defined sequence that may carry trailing gaps.
    """
    real = gapless_seq(seq)
    if len(real) < 3 or real[-3:] not in STOP_CODONS:
        return seq
    non_gap_idx = [i for i, c in enumerate(seq) if c != "-"]
    remove = set(non_gap_idx[-3:])
    if preserve_width:
        return "".join("-" if i in remove else c for i, c in enumerate(seq))
    return "".join(c for i, c in enumerate(seq) if i not in remove)


def orf_check(seq):
    # ORF status is assessed on the biological sequence, not on alignment gaps.
    # N is kept in the coordinate system and length count; exact STOP codons are
    # detected only when all three bases are known.
    real = gapless_seq(seq)
    n = len(real)
    starts = n >= 3 and real[:3] == "ATG"
    mod3 = n % 3 == 0
    terminal = real[-3:] if n >= 3 and mod3 else "NA"
    has_term = terminal in STOP_CODONS
    premature = []
    full = n // 3

    # Important: absence of a terminal STOP is not treated as ORF failure.
    # Many curated CDS FASTA files intentionally omit the final STOP codon
    # because downstream alignment/ASR programs can mishandle terminal '*'.
    # Therefore, a sequence is coding-complete if it starts with ATG, has a
    # length divisible by three, and lacks internal in-frame STOP codons.
    # If a terminal STOP is present, it is excluded from the internal-stop scan.
    # If no terminal STOP is present, all full codons are scanned.
    scan_n = full - 1 if (mod3 and has_term) else full
    for i in range(scan_n):
        codon = real[i*3:i*3+3]
        if codon in STOP_CODONS:
            premature.append((i+1, i*3+1, i*3+3, codon))
    complete = starts and mod3 and len(premature) == 0
    return {
        "sequence_length": n,
        "starts_with_atg": starts,
        "length_multiple_of_3": mod3,
        "terminal_codon": terminal,
        "has_terminal_stop": has_term,
        "premature_stop_count": len(premature),
        "premature_stop_codons": ",".join(str(x[0]) for x in premature) if premature else "NA",
        "complete_orf": complete,
        "premature_details": premature,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--alignment-mode", choices=["perform", "defined"], default="perform",
                    help="defined: the FASTA is an authoritative codon alignment; the terminal "
                         "stop is removed by gapping in place so every column is preserved.")
    args = ap.parse_args()

    gene = args.gene
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    preserve_width = args.alignment_mode == "defined"
    recs = {r.id.split()[0]: clean_seq(str(r.seq)) for r in SeqIO.parse(args.fasta, "fasta")}
    recs = {name: strip_terminal_stop(seq, preserve_width=preserve_width) for name, seq in recs.items()}
    tree = Phylo.read(args.tree, "newick")
    tree_tips = {t.name for t in tree.get_terminals()}
    fasta_names = set(recs)
    common_set = fasta_names & tree_tips

    if not common_set:
        raise SystemExit("No matching names between FASTA IDs and tree tip labels.")

    dropped = []
    for s in sorted(fasta_names - tree_tips):
        dropped.append({"gene": gene, "species": s, "reason": "fasta_only"})
    for s in sorted(tree_tips - fasta_names):
        dropped.append({"gene": gene, "species": s, "reason": "tree_only"})

    write_tsv(dropped, outdir / f"00_{gene}.dropped_species.tsv", ["gene", "species", "reason"])

    for tip in list(tree.get_terminals()):
        if tip.name not in common_set:
            tree.prune(tip)
    Phylo.write(tree, outdir / f"00_{gene}.common_species.tree", "newick")

    # Preserve the tip order of the rooted, pruned user tree, matching how
    # 05_plot_events.R (ape::plot.phylo, direction="rightwards") actually
    # renders it: the FIRST tip in plain left-to-right preorder ends up drawn
    # at the BOTTOM of the figure and the LAST at the TOP, so writing files in
    # that plain order reads bottom-to-top next to the rendered tree --
    # backwards from how a person reads a file top to bottom. Visiting each
    # node's own children right-to-left instead (still parent-before-children,
    # so still a real preorder) produces the order that reads top-to-bottom
    # in the same order as the plot. Pensieve never alphabetically reorders
    # biological sequences for MSA output either way.
    def terminals_top_to_bottom(clade):
        if not clade.clades:
            return [clade]
        tips = []
        for child in reversed(clade.clades):
            tips.extend(terminals_top_to_bottom(child))
        return tips
    common = [tip.name for tip in terminals_top_to_bottom(tree.root)]

    with open(outdir / f"00_{gene}.common_species.fasta", "w") as out:
        for s in common:
            out.write(f">{s}\n")
            seq = recs[s]
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

    # MACSE receives gapless reconstructed CDS even when the user supplied a
    # curated alignment for --alignment defined.  This keeps diagnostic MACSE
    # behaviour independent from manual alignment gaps.
    with open(outdir / f"00_{gene}.common_species.gapless_for_macse.fasta", "w") as out:
        for s in common:
            out.write(f">{s}\n")
            seq = gapless_seq(recs[s])
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

    write_tsv([
        {"gene": gene, "phylogenetic_tip_rank": i, "species": s}
        for i, s in enumerate(common, start=1)
    ], outdir / f"00_{gene}.phylogenetic_tip_order.tsv",
       ["gene", "phylogenetic_tip_rank", "species"])

    status_rows = []
    failure_rows = []
    complete_rows = []
    orf_by_species = {}

    for s in common:
        ch = orf_check(recs[s])
        orf_by_species[s] = ch
        row = {"gene": gene, "species": s}
        for k in ["sequence_length", "starts_with_atg", "length_multiple_of_3", "terminal_codon", "has_terminal_stop", "premature_stop_count", "premature_stop_codons", "complete_orf"]:
            row[k] = ch[k]
        status_rows.append(row)
        if ch["complete_orf"]:
            complete_rows.append({"gene": gene, "species": s, "sequence_length": ch["sequence_length"]})
        if not ch["starts_with_atg"]:
            failure_rows.append({"gene": gene, "species": s, "failure_type": "missing_start_ATG", "codon_position": 1, "nt_start": 1, "nt_end": 3, "codon": gapless_seq(recs[s])[:3], "details": "gapless_sequence"})
        if not ch["length_multiple_of_3"]:
            failure_rows.append({"gene": gene, "species": s, "failure_type": "length_not_multiple_of_3", "codon_position": "NA", "nt_start": "NA", "nt_end": "NA", "codon": "NA", "details": f"length={ch['sequence_length']}"})
        # Missing terminal STOP is deliberately not a failure.  The terminal
        # STOP is often absent in collaborator-provided CDS files.
        for codpos, nt1, nt2, codon in ch["premature_details"]:
            failure_rows.append({"gene": gene, "species": s, "failure_type": "premature_in_frame_stop", "codon_position": codpos, "nt_start": nt1, "nt_end": nt2, "codon": codon, "details": "raw_unaligned_sequence"})

    # Split the gapless CDS by ORF completeness so MACSE can be told which
    # sequences it is allowed to trust. Complete ORFs go to `-seq` (reliable);
    # everything else goes to `-seq_lr` (less reliable, MACSE's own documented
    # slot for pseudogenes). A frameshift proposed inside a reliable sequence
    # then costs `-fs`, which step 01 sets prohibitively high, so MACSE cannot
    # "fix" an intact reading frame in order to make a broken lineage align
    # more cheaply -- the failure mode the removed conserved-block heuristic
    # was patching after the fact.
    complete_species = [s for s in common if orf_by_species[s]["complete_orf"]]
    incomplete_species = [s for s in common if not orf_by_species[s]["complete_orf"]]

    def write_gapless(names, path):
        with open(path, "w") as out:
            for s in names:
                seq = gapless_seq(recs[s])
                out.write(f">{s}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i + 80] + "\n")

    write_gapless(complete_species, outdir / f"00_{gene}.complete_seqs.fa")
    write_gapless(incomplete_species, outdir / f"00_{gene}.incomplete_seqs.fa")
    print(f"ORF split for {gene}: {len(complete_species)} complete ORF sequence(s) -> "
          f"00_{gene}.complete_seqs.fa (MACSE -seq); {len(incomplete_species)} incomplete -> "
          f"00_{gene}.incomplete_seqs.fa (MACSE -seq_lr)")
    if not complete_species:
        print(f"[WARN] {gene}: no complete-ORF sequence; step 01 will pass every sequence to MACSE "
              f"as reliable (-seq) because -seq_lr alone is not a valid MACSE input.", file=sys.stderr)

    write_tsv(status_rows, outdir / f"00_{gene}.orf_status.tsv", ["gene", "species", "sequence_length", "starts_with_atg", "length_multiple_of_3", "terminal_codon", "has_terminal_stop", "premature_stop_count", "premature_stop_codons", "complete_orf"])
    write_tsv(failure_rows, outdir / f"00_{gene}.orf_failures.tsv", ["gene", "species", "failure_type", "codon_position", "nt_start", "nt_end", "codon", "details"])
    write_tsv(complete_rows, outdir / f"00_{gene}.complete_orf_species.tsv", ["gene", "species", "sequence_length"])

    print(f"Finished step00 for {gene}; common species: {len(common)}; "
          f"complete ORFs: {len(complete_species)}; incomplete: {len(incomplete_species)}; "
          f"coordinate system: canonical alignment")


if __name__ == "__main__":
    main()
