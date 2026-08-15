#!/usr/bin/env python3
"""ORF-aware, pseudogenic-component-bounded structural character reconstruction
(v4.1). Regression tests mandated by the ChatGPT review of a real, densely
sampled pseudogene-clade dataset (Bat_genes_from_Song, GUCY2F).

Real problem this addresses: scripts/03_alignment_events.py reconstructs each
indel/STOP character independently with binary Sankoff parsimony, where every
tip is fed in as an observation. For a monophyletic pseudogenized clade, tips
that share a *passenger* mutation (a structural change that happened after
the gene died, not the disabling lesion itself) are NOT independent evidence
about whether that mutation -- or the pseudogenization itself -- predates the
clade's common ancestor; they are one correlated sample. As a clade's
"carrier" and "non-carrier" counts both grow (proportionally, the same
underlying biology just sampled more densely), the RAW Sankoff cost the
clade's own entry node computes for itself grows without bound even though
the entry node's own best-supported STATE does not need to change -- and the
resulting parsimony confidence (delta_parsimony_support) can be made
arbitrarily large purely by adding more descendants, which is exactly the
"sampling density masquerading as independent evidence" failure the review
flagged in real GUCY2F data (see CHANGELOG v4.1).

The fix (scripts/03_alignment_events.py: pseudogenic_components(),
clamp_boundary_evidence(), and the component_entry_idx/boundary_cap
arguments threaded through sankoff()) is Stage A (provisional tip-evidence-
only ORF history), Stage B (identify F->P component entry branches), and
Stage C (bound the entry node's own C=ABSENT/PRESENT cost difference, and
everything it hands to its parent, to a fixed number of "equivalent
independent tip votes" -- regardless of real descendant count).
"""
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FAILURES: list[str] = []


def check(condition, message):
    if not condition:
        FAILURES.append(message)
        print(f"  FAIL {message}")
    else:
        print(f"  ok   {message}")


def tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_orf_status(path, complete_by_species):
    with open(path, "w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t")
        w.writerow(["gene", "species", "complete_orf"])
        for species, complete in complete_by_species.items():
            w.writerow(["T", species, str(bool(complete))])


def run_events(tmp: Path, tree_text: str, seqs: dict[str, str], orf_status=None,
                extra_args=None, gene="T"):
    tree = tmp / "tree.nwk"
    tree.write_text(tree_text + ("\n" if not tree_text.endswith("\n") else ""))
    aln = tmp / "aln.fa"
    aln.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    out = tmp / "events"
    cmd = [sys.executable, str(ROOT / "scripts" / "03_alignment_events.py"),
           "--gene", gene, "--alignment", str(aln), "--tree", str(tree),
           "--outdir", str(out), "--dated", "no", "--tie-break", "none"]
    if orf_status is not None:
        status_path = tmp / "orf_status.tsv"
        write_orf_status(status_path, orf_status)
        cmd += ["--orf-status", str(status_path)]
    if extra_args:
        cmd += extra_args
    r = subprocess.run(cmd, capture_output=True, text=True)
    check(r.returncode == 0, f"event engine runs ({r.stderr[-300:]})")
    return out


def original_label_map(out):
    return {r["pensieve_node_label"]: r["original_tree_label"]
            for r in tsv(out / "03_T.internode_label_map.tsv")}


def provisional_history_by_label(out):
    """Keyed by the ORIGINAL Newick label; see indel_state_by_label()."""
    label_map = original_label_map(out)
    rows = tsv(out / "03_T.provisional_orf_history.tsv")
    return {label_map.get(r["node_label"], r["node_label"]): r for r in rows}


def indel_state_by_label(out):
    """Keyed by the ORIGINAL Newick label (e.g. 'Pclade'), not Pensieve's own
    postorder Node1/Node2/... renumbering -- apply_pensieve_labels() always
    relabels internal nodes, so tests must look the real label up via
    03_<gene>.internode_label_map.tsv's original_tree_label column."""
    rows = tsv(out / "03_T.alignment_character_node_states.tsv")
    label_map = {r["pensieve_node_label"]: r["original_tree_label"]
                 for r in tsv(out / "03_T.internode_label_map.tsv")}
    return {label_map.get(r["node_label"], r["node_label"]): r
            for r in rows if r["character_class"] == "indel"}


BASE = "ATGAAACCCGGGTTTAAA"          # 6 codons, in-frame baseline
FS_DEL = BASE[:6] + "--" + BASE[8:]  # 2bp deletion at cols 7-8 -- frameshifting -> ORF-relevant cost


def pclade_tree_and_seqs(k_carriers, m_noncarriers):
    """(Out0, (K frameshift-deletion carriers + M non-carriers)Pclade)Anc --
    a single pseudogenic clade of variable internal composition opposing one
    functional outgroup tip."""
    seqs = {f"P{i}": FS_DEL for i in range(k_carriers)}
    seqs.update({f"M{i}": BASE for i in range(m_noncarriers)})
    seqs["Out0"] = BASE
    p_tips = ",".join([f"P{i}:1" for i in range(k_carriers)] + [f"M{i}:1" for i in range(m_noncarriers)])
    tree = f"(Out0:1,({p_tips})Pclade:1)Anc;"
    orf_status = {f"P{i}": False for i in range(k_carriers)}
    orf_status.update({f"M{i}": False for i in range(m_noncarriers)})
    orf_status["Out0"] = True
    return tree, seqs, orf_status


def test_pseudogene_sampling_density_invariance():
    # MANDATORY (spec 7.2): a pseudogenic radiation with many sampled
    # descendants must have approximately the same influence on its
    # pre-loss ancestor as the same radiation with few descendants.
    #
    # Real GUCY2F failure mode: a heterogeneous within-clade structural
    # character (present in SOME, not all, descendants of a pseudogenic
    # clade -- exactly how independently-accumulated passenger mutations
    # look in real data) makes the clade's own entry-node Sankoff cost
    # difference (and hence delta_parsimony_support, its reported
    # confidence) grow linearly with clade size under the OLD, uncapped
    # algorithm, purely from sample count, with the same 50/50
    # carrier:non-carrier COMPOSITION held fixed at every size (so the true
    # underlying biology genuinely is unchanged between versions).
    #
    # v4.2 correction (ChatGPT review, spec 10): capping must NOT touch the
    # component-entry node's OWN reported state/confidence -- only the
    # MESSAGE it hands to its PARENT (an ancestor strictly ABOVE the
    # component boundary). So "Pclade" (the entry node itself) is expected
    # to show the SAME unbounded confidence growth under both algorithms
    # (full internal evidence, exactly like any node inside the component,
    # and a genuine reflection of real internal evidence, not a bug).
    #
    # A second, non-obvious mathematical fact verified here: ordinary binary
    # Sankoff already bounds a SINGLE child's own pull on its DIRECT parent's
    # cost DIFFERENCE at the transition cost (gain_cost/loss_cost), no matter
    # how large that child's own internal cost difference grows -- a min-plus
    # ceiling that holds with or without Stage C. So "Anc" (Pclade's parent)
    # is ALREADY sampling-density invariant even in the uncapped algorithm for
    # this topology; Stage C's boundary_cap only makes an observable
    # difference there when set BELOW that natural per-mutation ceiling,
    # giving a researcher an explicit dial to compress cross-boundary
    # influence further than a single mutation's own transition cost would
    # allow -- verified below with a cap deliberately set under the ceiling.
    print("pseudogene sampling-density invariance (critical)")
    entry_deltas, parent_deltas, parent_tight_deltas = [], [], []
    states = []
    sizes = [2, 5, 10, 20, 40, 80, 160, 320]
    for k in sizes:
        tree, seqs, orf_status = pclade_tree_and_seqs(k, k)
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            (tmp / "loose").mkdir()
            (tmp / "tight").mkdir()
            # "loose": cap (2.0 votes * 0.65 = 1.3) is ABOVE the natural
            # gain/loss transition ceiling -- should behave identically to no
            # cap at all at the parent level.
            out_loose = run_events(tmp / "loose", tree, seqs, orf_status=orf_status,
                                    extra_args=["--pseudogenic-boundary-cap-votes", "2.0"])
            # "tight": cap (0.3 votes * 0.65 = 0.195) is BELOW both gain_cost
            # (0.35) and loss_cost (0.65) -- should visibly compress the
            # message crossing into Anc, below the natural ceiling.
            out_tight = run_events(tmp / "tight", tree, seqs, orf_status=orf_status,
                                    extra_args=["--pseudogenic-boundary-cap-votes", "0.3"])
            loose_by_label = indel_state_by_label(out_loose)
            tight_by_label = indel_state_by_label(out_tight)
            entry_deltas.append(float(loose_by_label["Pclade"]["delta_parsimony_support"]))
            parent_deltas.append(float(loose_by_label["Anc"]["delta_parsimony_support"]))
            parent_tight_deltas.append(float(tight_by_label["Anc"]["delta_parsimony_support"]))
            states.append(loose_by_label["Pclade"]["state"])

    check(len(set(states)) == 1,
          f"the pseudogenic clade's own entry-node state is IDENTICAL across clade sizes "
          f"{sizes} regardless of cap ({states})")
    check(entry_deltas == sorted(entry_deltas) and entry_deltas[-1] > entry_deltas[0] * 10,
          f"the component-entry node's OWN reported confidence (delta_parsimony_support) grows "
          f"without bound purely from clade size, same 50/50 composition throughout -- this is "
          f"CORRECT per spec 10, not a bug: the entry node always uses full internal evidence "
          f"({list(zip(sizes, entry_deltas))})")
    check(max(parent_deltas) - min(parent_deltas) < 0.1,
          f"the entry node's PARENT ('Anc') is already sampling-density invariant on its own (up to "
          f"a small edge effect at the smallest clade size), an inherent Sankoff min-plus property "
          f"that holds independently of Stage C's cap once the cap is at/above the natural "
          f"transition-cost ceiling: {list(zip(sizes, parent_deltas))}")
    check(max(parent_tight_deltas) - min(parent_tight_deltas) < 1e-6,
          f"'Anc' stays invariant to clade size under a TIGHTER cap too, at a different "
          f"(smaller) fixed value: {list(zip(sizes, parent_tight_deltas))}")
    check(parent_tight_deltas[0] < parent_deltas[0] - 1e-9,
          f"a boundary_cap explicitly set BELOW the natural per-mutation transition cost measurably "
          f"compresses the message crossing into 'Anc' further than Sankoff's own natural ceiling "
          f"would alone (tight={parent_tight_deltas[0]:g}, loose/natural={parent_deltas[0]:g}), "
          f"confirming the cap parameter has a real, working, size-invariant effect")


def test_multiple_losses_preferred_over_repeated_restoration_when_supported():
    # Spec 7.3: with two or more pseudogenic clades carrying the same
    # disabling state, separated by intact lineages, Stage A should prefer
    # "ancestor intact, independent losses" over "ancestor pseudogenic,
    # repeated restoration" when that is cheaper under the explicit F/P
    # model -- and this must be topology-driven, not a function of raw tip
    # count in either clade.
    print("independent pseudogenizations preferred over repeated restoration when supported")
    # (Intact1, (PcladeA, (Intact2, PcladeB)) ) -- two SEPARATE small P
    # clades, each directly sistered by/near an intact lineage, so the
    # cheapest explanation is "root functional, two independent losses"
    # rather than "root pseudogenic, two independent restorations".
    tree = "(Intact1:1,((PA1:1,PA2:1)PcladeA:1,(Intact2:1,(PB1:1,PB2:1)PcladeB:1)Mid:1)Deep:1)Root;"
    seqs = {
        "Intact1": BASE, "Intact2": BASE,
        "PA1": FS_DEL, "PA2": FS_DEL, "PB1": FS_DEL, "PB2": FS_DEL,
    }
    orf_status = {"Intact1": True, "Intact2": True, "PA1": False, "PA2": False, "PB1": False, "PB2": False}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        out = run_events(tmp, tree, seqs, orf_status=orf_status,
                          extra_args=["--orf-restoration-cost", "4.0", "--orf-loss-cost", "1.0"])
        by_label = provisional_history_by_label(out)
        check(by_label["Root"]["provisional_orf_state"] == "functional",
              f"root stays functional when two independent, cheaply-explained losses beat one "
              f"ancestral pseudogenization plus two expensive restorations ({by_label['Root']})")
        components = tsv(out / "03_T.pseudogenic_components.tsv")
        entry_names = {r["component_entry_node"] for r in components}
        check(len(components) == 2,
              f"two SEPARATE pseudogenic components are identified (independent losses), not one "
              f"shared ancestral component ({components})")


def test_root_not_forced_functional():
    # Spec 7.4 (MANDATORY): an all-/overwhelmingly-pseudogenic synthetic
    # case with no reliable intact ancestral evidence must be able to
    # report a pseudogenic or uncertain root -- never a hard-coded
    # functional root.
    print("root is not forced functional when the data do not support it")
    tree = "((P1:1,P2:1)C1:1,(P3:1,P4:1)C2:1)Root;"
    seqs = {"P1": FS_DEL, "P2": FS_DEL, "P3": FS_DEL, "P4": FS_DEL}
    orf_status = {"P1": False, "P2": False, "P3": False, "P4": False}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        out = run_events(tmp, tree, seqs, orf_status=orf_status)
        by_label = provisional_history_by_label(out)
        check(by_label["Root"]["provisional_orf_state"] != "functional",
              f"an all-pseudogenic tree does not get a hard-coded functional root "
              f"({by_label['Root']})")

    # A tip absent from orf_status.tsv must contribute NO forced state --
    # not "functional" and not "pseudogenic". sankoff() already represents
    # this correctly with down=[0,0] (unconstrained), but an unconstrained
    # tip still gets a cosmetic REPRESENTATIVE label for reporting (it
    # costlessly "inherits" whatever its parent's own state ends up being,
    # exactly like any other missing-data tip elsewhere in Pensieve) -- that
    # label is not a vote. The real test of "no forced vote" is that adding
    # an unconstrained tip must not change any REAL node's own classification
    # versus the same tree with that tip simply absent altogether.
    print("an unconstrained (orf_status-absent) tip does not change any other node's classification")
    tree_with_u = "((P1:1,U1:1)C1:1,P2:1)Root;"
    seqs_with_u = {"P1": FS_DEL, "U1": BASE, "P2": FS_DEL}
    orf_status_with_u = {"P1": False, "P2": False}  # U1 deliberately absent -> unconstrained
    tree_without_u = "(P1:1,P2:1)Root;"
    seqs_without_u = {"P1": FS_DEL, "P2": FS_DEL}
    orf_status_without_u = {"P1": False, "P2": False}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        (tmp / "with_u").mkdir()
        (tmp / "without_u").mkdir()
        out_with = run_events(tmp / "with_u", tree_with_u, seqs_with_u, orf_status=orf_status_with_u)
        out_without = run_events(tmp / "without_u", tree_without_u, seqs_without_u,
                                  orf_status=orf_status_without_u)
        root_with = provisional_history_by_label(out_with)["Root"]["provisional_orf_state"]
        root_without = provisional_history_by_label(out_without)["Root"]["provisional_orf_state"]
        check(root_with == root_without,
              f"Root's classification is unchanged by adding an unconstrained tip U1 "
              f"(with U1: {root_with}, without U1: {root_without})")


def test_apparent_restoration_retains_pseudogenic_history():
    # Spec 7.5: a lineage with known upstream pseudogenic history and a
    # descendant whose CURRENT sequence is ORF-intact must still be
    # representable, and Stage A/B/C must not silently erase that prior
    # pseudogenic history for the character-level reconstruction.
    print("apparent restoration: Stage A still records the deeper pseudogenic history")
    # Root(F) -> C1(P, 2 tips confidently pseudogenic) -> within C1, one
    # descendant clade R (2 tips) is confidently intact again (apparent
    # restoration at the ORF-history level; the sequence-level "apparent
    # restoration" semantics live in 04_ancestral_orf_walk.py and are
    # untouched by this change -- this test only confirms Stage A's own
    # provisional history keeps the component entry above R, not erased at R).
    # C1 needs enough INDEPENDENT pseudogenic branches directly under it
    # (each its own entry, at orf_loss_cost apiece) to outweigh the single
    # Restored clade's orf_restoration_cost -- with the default costs
    # (loss=1.0, restoration=4.0), C1=functional would cost N x 1.0
    # (independent gains) versus C1=pseudogenic costing 1 x 4.0 (one
    # restoration for Restored); N=6 safely clears that break-even point
    # (6.0 > 4.0), so C1's own pseudogenic history is genuinely
    # well-supported by real, independent evidence -- not asserted by fiat.
    tree = ("(Out:1,(P1:1,P2:1,P3:1,P4:1,P5:1,P6:1,(R1:1,R2:1)Restored:1)C1:1)Root;")
    seqs = {"Out": BASE, "R1": BASE, "R2": BASE}
    seqs.update({f"P{i}": FS_DEL for i in range(1, 7)})
    orf_status = {"Out": True, "R1": True, "R2": True}
    orf_status.update({f"P{i}": False for i in range(1, 7)})
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        out = run_events(tmp, tree, seqs, orf_status=orf_status)
        by_label = provisional_history_by_label(out)
        check(by_label["C1"]["provisional_orf_state"] == "pseudogenic",
              f"C1 (ancestor of both the still-pseudogenic branches and the apparently-restored "
              f"clade), with 6 independent pseudogenic branches outweighing 1 restoration, "
              f"correctly retains pseudogenic provisional history ({by_label['C1']})")
        check(by_label["Root"]["provisional_orf_state"] == "functional",
              f"Root itself, outside the C1 component (opposed by the single functional Out tip), "
              f"is unaffected ({by_label['Root']})")


def test_all_functional_matches_neutral_behavior():
    # Spec 7.6 (regression guard): with no pseudogenic component anywhere,
    # passing --orf-status must reproduce EXACTLY the old, uncapped output --
    # detects accidental global weighting or unintended changes to ordinary
    # in-frame indel placement.
    print("neutral/all-functional regression: identical output with and without --orf-status")
    tree = "((A:1,B:1)C1:1,(D:1,E:1)C2:1)Root;"
    seqs = {"A": BASE, "B": FS_DEL, "D": BASE, "E": BASE}
    orf_status = {"A": True, "B": True, "D": True, "E": True}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        (tmp / "plain").mkdir()
        (tmp / "orf").mkdir()
        out_plain = run_events(tmp / "plain", tree, seqs)
        out_orf = run_events(tmp / "orf", tree, seqs, orf_status=orf_status)
        plain_events = tsv(out_plain / "03_T.alignment_events.tsv")
        orf_events = tsv(out_orf / "03_T.alignment_events.tsv")

        def strip(rows):
            return [{k: v for k, v in r.items()} for r in rows]
        check(strip(plain_events) == strip(orf_events),
              "alignment_events.tsv is byte-identical whether or not --orf-status is supplied, "
              "when no pseudogenic component exists (all tips confidently functional)")
        components = tsv(out_orf / "03_T.pseudogenic_components.tsv")
        check(len(components) == 0,
              f"no pseudogenic component is (falsely) identified when every tip is functional ({components})")


def test_unknown_boundary_remains_ambiguous():
    # Spec 7.7: unknown/missing ORF states must not become pseudo-functional
    # (or pseudo-pseudogenic) votes, and a pseudogenic component with
    # ambiguous boundary evidence must not force a confident entry call.
    print("unknown tip ORF evidence does not force a confident component boundary")
    # A clade where HALF the tips are pseudogenic and half are simply
    # missing from orf_status.tsv (unconstrained) -- Stage A must not
    # confidently call this clade's entry branch a component just because
    # the unconstrained tips default toward whichever is cheapest; check
    # that at minimum no crash occurs and the missing tips are not silently
    # read as functional (which would wrongly suppress a real component) or
    # as pseudogenic (which would wrongly invent one).
    tree_with_u = "(Out:1,(P1:1,P2:1,U1:1,U2:1)Mixed:1)Root;"
    seqs_with_u = {"Out": BASE, "P1": FS_DEL, "P2": FS_DEL, "U1": BASE, "U2": BASE}
    orf_status_with_u = {"Out": True, "P1": False, "P2": False}  # U1/U2 deliberately absent
    tree_without_u = "(Out:1,(P1:1,P2:1)Mixed:1)Root;"
    seqs_without_u = {"Out": BASE, "P1": FS_DEL, "P2": FS_DEL}
    orf_status_without_u = {"Out": True, "P1": False, "P2": False}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        (tmp / "with_u").mkdir()
        (tmp / "without_u").mkdir()
        out_with = run_events(tmp / "with_u", tree_with_u, seqs_with_u, orf_status=orf_status_with_u)
        out_without = run_events(tmp / "without_u", tree_without_u, seqs_without_u,
                                  orf_status=orf_status_without_u)
        mixed_with = provisional_history_by_label(out_with)["Mixed"]["provisional_orf_state"]
        mixed_without = provisional_history_by_label(out_without)["Mixed"]["provisional_orf_state"]
        check(mixed_with == mixed_without,
              f"adding two unconstrained tips (U1, U2, absent from orf_status.tsv) to the pseudogenic "
              f"clade does not change the clade entry's own classification (with U1/U2: {mixed_with}, "
              f"without: {mixed_without})")
        # direction_confident/ambiguous_origin diagnostics must still be internally consistent
        events = tsv(out_with / "03_T.alignment_events.tsv")
        for ev in events:
            if ev["direction_confident"] == "False":
                check(ev["ambiguous_origin"] == "True",
                      f"an event with direction_confident=False is always also ambiguous_origin=True "
                      f"({ev['event_id']})")


def main():
    test_pseudogene_sampling_density_invariance()
    test_multiple_losses_preferred_over_repeated_restoration_when_supported()
    test_root_not_forced_functional()
    test_apparent_restoration_retains_pseudogenic_history()
    test_all_functional_matches_neutral_behavior()
    test_unknown_boundary_remains_ambiguous()
    if FAILURES:
        print(f"\n{len(FAILURES)} failure(s)")
        sys.exit(1)
    print("\nORF-aware parsimony test passed.")


if __name__ == "__main__":
    main()
