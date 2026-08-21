#!/usr/bin/env python3
"""Functional-shared ancestral indel pre-pass (Pensieve stage `diagnostics`).

Runs immediately after the alignment exists and before anything is
reconstructed:

  * `--alignment perform`  -- after `01_run_macse_and_extract_events.py`, on
    `01_<gene>.macse_NT.fasta` (MACSE's `!` partial-codon placeholder rendered
    as `-`, exactly as step 02 will render the canonical native view, so the
    columns here and downstream are the same columns).
  * `--alignment defined`  -- after `00_prune_and_check_orf.py`, on
    `00_<gene>.common_species.fasta`, which step 02 uses unchanged.

It codes the tip alignment into breakpoint-defined indel characters with the
SAME functions the event reconstruction uses (`cluster_indel_characters`,
`decompose_run`, `indel_tip_states` in `03_alignment_events.py`, imported here
rather than copied), then compares, character by character, which
complete-ORF (functional) lineages carry the indel against which pseudogenized
lineages do.

The biology this encodes: functional lineages are trustworthy witnesses of the
ancestral functional sequence -- several INDEPENDENT functional lineages
acquiring exactly the same indel convergently is very close to impossible, so
an indel the functional lineages share was already present in their common
functional ancestor. Pseudogenized lineages are not witnesses in that sense: a
dead gene accumulates arbitrary indels, so any number of pseudogenized lineages
can independently gain or destroy any given indel, and their agreement is not
evidence of ancestry.

The result is written to `01b_<gene>.functional_shared_indels.tsv` -- every
indel character with the evidence and the verdict -- and handed to
`03_alignment_events.py --functional-shared-indels`, which fixes the flagged
indels as ancestral when it reconstructs the whole alignment (tips plus the
PAML ancestral sequences integrated downstream).
"""
import argparse
import importlib.util
from pathlib import Path

_EVENTS = Path(__file__).with_name("03_alignment_events.py")
_spec = importlib.util.spec_from_file_location("pensieve_alignment_events", _EVENTS)
ae = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ae)


def main():
    parser = argparse.ArgumentParser(
        description="Identify indels shared by the functional lineages (ancestral) before reconstruction.")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--alignment", required=True,
                        help="Diagnostics-stage alignment: 01_<gene>.macse_NT.fasta (perform) or "
                             "00_<gene>.common_species.fasta (defined).")
    parser.add_argument("--tree", required=True, help="00_<gene>.common_species.tree")
    parser.add_argument("--orf-status", required=True,
                        help="00_<gene>.orf_status.tsv -- the tip-level complete_orf calls that define "
                             "which lineages are functional.")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--breakpoint-tolerance", type=int, default=0,
                        help="Must match the value used by 03_alignment_events.py so both passes code "
                             "the same characters.")
    parser.add_argument("--min-functional-witnesses", type=int, default=2,
                        help="Minimum number of independent functional carriers; 0 disables the rule.")
    args = parser.parse_args()

    gene = args.gene
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    alignment, _order, aln_len = ae.read_fasta(args.alignment)
    # MACSE's '!' placeholder is not a residue and not an indel of its own; step
    # 02 renders it '-' in the canonical native view, so render it the same way
    # here or the gap runs (and therefore the characters) would differ.
    alignment = {name: seq.replace("!", "-") for name, seq in alignment.items()}

    root = ae.parse_newick(args.tree)
    root, _collapsed = ae.collapse_unifurcations(root)
    root, _labels = ae.apply_pensieve_labels(root)
    nodes_post = list(ae.iter_nodes(root, "postorder"))
    tips = [n.name for n in nodes_post if n.is_tip]

    missing = sorted(set(tips) - set(alignment))
    if missing:
        raise SystemExit(f"{len(missing)} tree tip(s) absent from the alignment: {missing[:10]}")

    orf_rows = ae.read_tsv(args.orf_status)
    tip_orf_states = ae.read_tip_orf_states(orf_rows, tips)
    functional_tips = frozenset(t for t in tips if tip_orf_states.get(t) == ae.ABSENT)
    tips_under = ae.subtree_tip_sets(nodes_post)
    node_by_name = {n.name: n for n in nodes_post}

    print(f"[info] {gene}: {len(tips)} tips, alignment {aln_len} columns, "
          f"{len(functional_tips)} complete-ORF (functional) lineage(s), "
          f"{len(tips) - len(functional_tips)} pseudogenized/incomplete")
    if not functional_tips:
        print("[warn] no complete-ORF lineage: no indel can be witnessed as functional-shared")

    clusters, segments_by_tip, raw_runs = ae.cluster_indel_characters(
        alignment, tips, args.breakpoint_tolerance)

    rows = []
    for i, cluster in enumerate(clusters, start=1):
        start, end = cluster["start"], cluster["end"]
        states = ae.indel_tip_states(alignment, tips, segments_by_tip, raw_runs,
                                     start, end, args.breakpoint_tolerance)
        verdict = ae.functional_shared_indel_verdict(
            states, functional_tips, root, tips_under, node_by_name,
            args.min_functional_witnesses)
        length = end - start + 1
        rows.append({
            "gene": gene,
            "character_id": f"IND{i:04d}",
            "alignment_start": start,
            "alignment_end": end,
            "character_length": length,
            "length_mod_3": length % 3,
            "frame_effect": "in_frame" if length % 3 == 0 else "frameshift",
            "n_functional_present": verdict["n_functional_present"],
            "n_functional_absent": verdict["n_functional_absent"],
            "n_functional_unknown": verdict["n_functional_unknown"],
            "n_pseudogenized_present": verdict["n_pseudogenized_present"],
            "n_pseudogenized_absent": verdict["n_pseudogenized_absent"],
            "n_independent_functional_lineages": verdict["n_independent_functional_lineages"],
            "mrca_of_functional_carriers": verdict["mrca_of_functional_carriers"],
            "ancestral_functional_indel": verdict["is_ancestral"],
            "reason": verdict["reason"],
            "functional_carriers": verdict["functional_carriers"],
        })

    path = out / f"01b_{gene}.functional_shared_indels.tsv"
    ae.write_tsv(rows, path, ae.FUNCTIONAL_SHARED_HEADER)

    ancestral = [r for r in rows if r["ancestral_functional_indel"]]
    print(f"[result] {len(rows)} indel character(s) scored; {len(ancestral)} are shared by the "
          f"functional lineages and will be fixed as ancestral "
          f"({sum(1 for r in ancestral if r['frame_effect'] == 'frameshift')} frameshift, "
          f"{sum(1 for r in ancestral if r['frame_effect'] == 'in_frame')} in-frame)")
    for r in sorted(ancestral, key=lambda r: -r["n_functional_present"])[:15]:
        print(f"   {r['character_id']} {r['alignment_start']}-{r['alignment_end']} "
              f"({r['character_length']} bp, {r['frame_effect']}): functional "
              f"{r['n_functional_present']} carry / {r['n_functional_absent']} lack, "
              f"pseudogenized {r['n_pseudogenized_present']} carry / "
              f"{r['n_pseudogenized_absent']} lack, MRCA {r['mrca_of_functional_carriers']}")
    print(f"[done] {path}")


if __name__ == "__main__":
    main()
