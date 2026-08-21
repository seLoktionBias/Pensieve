#!/usr/bin/env python3
"""Functional-shared ancestral indel rule + event merging (v4.7).

Real problem this addresses (found by inspecting the real 103-species output
for the eight bat visual genes, CNGB3 and PDE6C most clearly): identical indel
events were reported scattered over many separate functional branches, so a
gene read as if it had been broken and resurrected again and again.

Two independent causes, both regression-tested here.

1. THE FUNCTIONAL-CONSENSUS GUARD. The v4.5/v4.6 rule fixed an indel as
   ancestral as soon as >=2 independent complete-ORF lineages carried it, with
   no regard for how many functional lineages definitely LACKED it. On real
   data that inverted the very implausibility the rule is built on: real PDE6C
   cols 107-108 had 7 functional carriers against 83 functional non-carriers
   and was pinned ancestral, producing 1 gain + 20 losses -- i.e. asserting 20
   independent, exact restorations of a 2 bp gap. Real CNGB3 3319-3324 and
   3355-3357 behaved the same way (4-8 carriers vs 20-24 non-carriers, 11-15
   loss events each). The rule now also requires the indel to be the functional
   CONSENSUS. With that guard the v4.6 frameshift-only restriction is no longer
   needed (it existed only to contain the un-guarded rule on real CNGA3
   409/472, which the consensus guard rejects directly), so a genuine
   functional-consensus in-frame indel -- real CNGB3 781-783, 16 functional
   carriers vs 11 -- is finally reconstructed as ancestral.

2. EVENT MERGING. Fragments of one real indel on one branch were left
   unmerged whenever (a) parsimony ambiguity made one fragment print as
   `ambiguous_indel_change` instead of `deletion` while its neighbour printed
   `deletion` (116 such contiguous same-direction pairs on real CNGB3 alone,
   ~280 across the eight genes), or (b) the fragments were separated by columns
   that were already gap in BOTH the parent and the child, i.e. inherited gap
   carrying no real sequence (real CNGB3 Node92->Dasypterus_ega, 2260-2289 +
   2293-2304 around an inherited 2290-2292 gap: one 42 bp deletion reported as
   two).
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


def write_inputs(tmp, tree_text, seqs, orf_status):
    tmp.mkdir(parents=True, exist_ok=True)
    tree = tmp / "tree.nwk"
    tree.write_text(tree_text + "\n")
    aln = tmp / "aln.fa"
    aln.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    status = tmp / "orf_status.tsv"
    write_orf_status(status, orf_status)
    return tree, aln, status


def run_prepass(tmp, tree, aln, status, gene="T"):
    out = tmp / "prepass"
    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "01b_functional_shared_indels.py"),
         "--gene", gene, "--alignment", str(aln), "--tree", str(tree),
         "--orf-status", str(status), "--outdir", str(out)],
        capture_output=True, text=True)
    check(r.returncode == 0, f"pre-pass runs ({r.stderr[-300:]})")
    return {(int(row["alignment_start"]), int(row["alignment_end"])): row
            for row in tsv(out / f"01b_{gene}.functional_shared_indels.tsv")}


def run_events(tmp, tree, aln, status, shared=None, gene="T", extra=None, out_name="events"):
    out = tmp / out_name
    cmd = [sys.executable, str(ROOT / "scripts" / "03_alignment_events.py"),
           "--gene", gene, "--alignment", str(aln), "--tree", str(tree),
           "--orf-status", str(status), "--outdir", str(out),
           "--dated", "no", "--tie-break", "none"]
    if shared is not None:
        cmd += ["--functional-shared-indels", str(shared)]
    if extra:
        cmd += extra
    r = subprocess.run(cmd, capture_output=True, text=True)
    check(r.returncode == 0, f"event engine runs ({r.stderr[-300:]})")
    return out


# 12 codons; a deletion is carved out of the middle so it is never terminal.
BASE = "ATG" + "AAACCCGGGTTTAAACCCGGGTTTAAACCCGGG" + "TAA"


def with_gap(start, length):
    """BASE with a `length` bp gap at 1-based column `start`."""
    return BASE[:start - 1] + "-" * length + BASE[start - 1 + length:]


def ladder(functional, pseudogenized):
    """A pectinate tree interleaving functional and pseudogenized tips so that
    any set of >=2 functional carriers is phylogenetically independent."""
    tips = []
    for i in range(max(len(functional), len(pseudogenized))):
        if i < len(functional):
            tips.append(functional[i])
        if i < len(pseudogenized):
            tips.append(pseudogenized[i])
    tree = tips[-1] + ":1"
    for name in reversed(tips[:-1]):
        tree = f"({name}:1,{tree}):1"
    return tree + ";"


def scenario(n_functional_carriers, n_functional_noncarriers, n_pseudo_carriers,
             n_pseudo_noncarriers, gap_start=13, gap_len=3):
    gapped, intact = with_gap(gap_start, gap_len), BASE
    seqs, orf, functional, pseudo = {}, {}, [], []
    for i in range(n_functional_carriers):
        seqs[f"Fc{i}"] = gapped; orf[f"Fc{i}"] = True; functional.append(f"Fc{i}")
    for i in range(n_functional_noncarriers):
        seqs[f"Fn{i}"] = intact; orf[f"Fn{i}"] = True; functional.append(f"Fn{i}")
    # Pseudogenized tips carry an unambiguous disabling lesion of their own
    # (an internal premature STOP) so that complete_orf is genuinely False.
    broken_gapped = gapped[:21] + "TAA" + gapped[24:]
    broken_intact = intact[:21] + "TAA" + intact[24:]
    for i in range(n_pseudo_carriers):
        seqs[f"Pc{i}"] = broken_gapped; orf[f"Pc{i}"] = False; pseudo.append(f"Pc{i}")
    for i in range(n_pseudo_noncarriers):
        seqs[f"Pn{i}"] = broken_intact; orf[f"Pn{i}"] = False; pseudo.append(f"Pn{i}")
    return ladder(functional, pseudo), seqs, orf


def verdict_for(tmp, n_fc, n_fn, n_pc, n_pn, gap_start=13, gap_len=3):
    tree, seqs, orf = scenario(n_fc, n_fn, n_pc, n_pn, gap_start, gap_len)
    t, a, s = write_inputs(tmp, tree, seqs, orf)
    rows = run_prepass(tmp, t, a, s)
    return rows.get((gap_start, gap_start + gap_len - 1))


print("functional consensus is required, not just >=2 functional carriers")
with tempfile.TemporaryDirectory() as d:
    tmp = Path(d)
    # Real PDE6C 107-108 shape: a few functional carriers, most functional
    # lineages definitely lack it.
    row = verdict_for(tmp / "a", 3, 12, 0, 4)
    check(row is not None, "the minority indel is scored as a character")
    if row:
        check(row["ancestral_functional_indel"] == "False",
              "an indel carried by a MINORITY of functional lineages (3 carry / 12 lack) is NOT ancestral")
        check(row["reason"] == "not_the_functional_consensus_state",
              f"the recorded reason names the consensus test (got {row['reason'] if row else 'NA'})")

with tempfile.TemporaryDirectory() as d:
    tmp = Path(d)
    row = verdict_for(tmp / "b", 12, 3, 0, 4)
    check(row is not None and row["ancestral_functional_indel"] == "True",
          "the same indel IS ancestral when the functional lineages share it (12 carry / 3 lack)")

print("pseudogenized carriers can never make a minority indel ancestral")
with tempfile.TemporaryDirectory() as d:
    tmp = Path(d)
    # Real CNGB3 3319/3355 shape: many pseudogenized carriers, functional
    # lineages mostly lack it. A dead gene accumulates arbitrary indels, so
    # their agreement is not evidence of ancestry.
    row = verdict_for(tmp / "c", 3, 12, 20, 2)
    check(row is not None and row["ancestral_functional_indel"] == "False",
          "20 pseudogenized carriers do not override 3-carry/12-lack among the functional lineages")
    if row:
        check(int(row["n_pseudogenized_present"]) >= 20,
              "the pseudogenized carrier count is still recorded for audit")

print("in-frame indels are eligible (the v4.6 frameshift-only restriction is gone)")
with tempfile.TemporaryDirectory() as d:
    tmp = Path(d)
    row = verdict_for(tmp / "d", 10, 3, 0, 4, gap_start=13, gap_len=3)
    check(row is not None and row["frame_effect"] == "in_frame",
          "the test character is in-frame")
    check(row is not None and row["ancestral_functional_indel"] == "True",
          "a functional-consensus IN-FRAME indel is ancestral (real CNGB3 781-783)")

print("a lineage-specific indel of one functional clade is not ancestral")
with tempfile.TemporaryDirectory() as d:
    tmp = Path(d)
    gapped = with_gap(13, 3)
    seqs = {"Fa": gapped, "Fb": gapped, "Fc": BASE, "Fd": BASE}
    orf = {k: True for k in seqs}
    tree = "((Fa:1,Fb:1):1,(Fc:1,Fd:1):1);"
    t, a, s = write_inputs(tmp, tree, seqs, orf)
    rows = run_prepass(tmp, t, a, s)
    row = rows.get((13, 15))
    check(row is not None and row["ancestral_functional_indel"] == "False",
          "two sister functional carriers forming a clean clade are one lineage-specific indel, not ancestral")

print("the ancestral verdict actually suppresses the dispersed per-lineage events")
with tempfile.TemporaryDirectory() as d:
    tmp = Path(d)
    tree, seqs, orf = scenario(8, 2, 0, 6)
    t, a, s = write_inputs(tmp, tree, seqs, orf)
    prepass_dir = tmp / "prepass"
    rows = run_prepass(tmp, t, a, s)
    shared = prepass_dir / "01b_T.functional_shared_indels.tsv"
    out_pinned = run_events(tmp, t, a, s, shared=shared, out_name="events_pinned")
    out_free = run_events(tmp, t, a, s, extra=["--min-functional-witnesses", "0"],
                          out_name="events_free")
    def gains(out):
        return sum(1 for r in tsv(out / "03_T.alignment_events.tsv")
                   if r["character_class"] == "indel"
                   and int(r["alignment_start"]) == 13
                   and r["event_type"] == "indel_gain")
    check(gains(out_pinned) <= 1,
          f"with the ancestral pin the shared indel is at most one gain (got {gains(out_pinned)})")
    check(gains(out_pinned) <= gains(out_free),
          f"the pin never increases the number of gain events (pinned={gains(out_pinned)}, free={gains(out_free)})")
    chars = {r["character_id"]: r for r in tsv(out_pinned / "03_T.alignment_characters.tsv")}
    check(any(r["functional_ancestral_indel"] == "True" for r in chars.values()),
          "the character table flags the indel as a functional ancestral indel")

print("pre-pass file and in-line recomputation agree")
with tempfile.TemporaryDirectory() as d:
    tmp = Path(d)
    tree, seqs, orf = scenario(8, 2, 3, 5)
    t, a, s = write_inputs(tmp, tree, seqs, orf)
    rows = run_prepass(tmp, t, a, s)
    shared = tmp / "prepass" / "01b_T.functional_shared_indels.tsv"
    with_file = run_events(tmp, t, a, s, shared=shared, out_name="events_file")
    file_events = (with_file / "03_T.alignment_events.tsv").read_text()
    inline = run_events(tmp, t, a, s, out_name="events_inline")
    inline_events = (inline / "03_T.alignment_events.tsv").read_text()
    check(file_events == inline_events,
          "reconstruction is byte-identical whether the verdicts come from the pre-pass file or are recomputed")

print("event merging: fragments of one branch's own indel are reported as one event")
with tempfile.TemporaryDirectory() as d:
    tmp = Path(d)
    # A carries a single contiguous 6 bp deletion at 13-18. B, on the far side
    # of the tree, independently carries only its first half -- which is what
    # splits A's one real deletion into two breakpoint characters (13-15 and
    # 16-18), both originating on A's own branch.
    seqs = {
        "A": BASE[:12] + "------" + BASE[18:],
        "C": BASE,
        "B": BASE[:12] + "---" + BASE[15:],
        "D": BASE,
    }
    orf = {k: True for k in seqs}
    t, a, s = write_inputs(tmp, "((A:1,C:1):1,(B:1,D:1):1);", seqs, orf)
    out = run_events(tmp, t, a, s)
    a_events = [r for r in tsv(out / "03_T.alignment_events.tsv")
                if r["character_class"] == "indel"
                and r["affected_tips"] == "A" and r["event_type"] == "indel_gain"]
    check(len(a_events) == 1,
          f"A's single contiguous 6 bp deletion is ONE event, not two fragments "
          f"(got {sorted((r['alignment_start'], r['alignment_end']) for r in a_events)})")
    if len(a_events) == 1:
        check(a_events[0]["event_length"] == "6" and a_events[0]["frame_effect"] == "in_frame",
              f"the merged event carries the true 6 bp in-frame length, not a spurious "
              f"3 bp fragment (got {a_events[0]['event_length']} bp, {a_events[0]['frame_effect']})")
        check("merged_contiguous_same_type" in a_events[0]["breakpoint_relationships"],
              "the merged event records the fragments it was built from")

print("event merging: the merge rules themselves (direct, on constructed event rows)")
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "pensieve_alignment_events", ROOT / "scripts" / "03_alignment_events.py")
ae = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ae)


def event_row(start, end, kind="indel_gain", interp="deletion", tips="A",
              branch="P->A", confident=True):
    length = end - start + 1
    return {
        "gene": "T", "event_id": f"T|IND{start:04d}|{kind[6:]}{start:04d}",
        "character_class": "indel", "event_type": kind,
        "biological_interpretation": interp,
        "alignment_start": start, "alignment_end": end, "event_length": length,
        "length_mod_3": length % 3,
        "frame_effect": "in_frame" if length % 3 == 0 else "frameshift",
        "origin_node": "A", "origin_is_tip": True, "parent_node": "P", "branch": branch,
        "shared_event": False, "n_affected_tips": len(tips.split(",")), "affected_tips": tips,
        "reversal_below_origin": False, "secondary_changes_below_origin": "NA",
        "root_state": "absent", "parsimony_score": "1", "delta_parsimony_support": "1",
        "ambiguous_origin": not confident, "direction_confident": confident,
        "parent_age": "NA", "child_age": "NA", "age_interval": "NA",
        "n_observed_present": 1, "n_observed_absent": 3, "n_unknown": 0,
        "observed_present_tips": tips, "breakpoint_relationships": "NA",
        "coordinate_system": "test",
    }


# (a) A confident fragment and an ambiguous fragment of the SAME direction are
#     one event -- the real CNGB3 Micronycteris_megalotis 1231-1236 /
#     1237-1242 / 1243-1248 shape, which the old label-based grouping left as
#     three separate events.
rows = ae.merge_contiguous_same_type_indel_events([
    event_row(1231, 1236),
    event_row(1237, 1242, interp="ambiguous_indel_change", confident=False),
    event_row(1243, 1248),
])
check(len(rows) == 1, f"confident+ambiguous+confident contiguous fragments merge into one event (got {len(rows)})")
if len(rows) == 1:
    check(int(rows[0]["alignment_start"]) == 1231 and int(rows[0]["alignment_end"]) == 1248
          and int(rows[0]["event_length"]) == 18,
          f"the merged span is the true 18 bp deletion 1231-1248 (got "
          f"{rows[0]['alignment_start']}-{rows[0]['alignment_end']}, {rows[0]['event_length']} bp)")
    check(rows[0]["biological_interpretation"] == "ambiguous_indel_change"
          and rows[0]["direction_confident"] is False,
          "merging never invents confidence: a cluster containing an ambiguous fragment stays ambiguous")

# (b) Opposite directions are never merged, however adjacent.
rows = ae.merge_contiguous_same_type_indel_events([
    event_row(100, 105),
    event_row(106, 110, kind="indel_loss", interp="insertion_or_restoration"),
])
check(len(rows) == 2, f"an adjacent deletion and insertion stay two separate events (got {len(rows)})")

# (c) Disjoint affected tips are never merged: two unrelated characters can
#     land on one branch by coincidence with no evidence they are one indel.
rows = ae.merge_contiguous_same_type_indel_events([
    event_row(100, 105, tips="A"),
    event_row(106, 110, tips="B"),
])
check(len(rows) == 2, f"contiguous fragments with disjoint affected tips stay separate (got {len(rows)})")

# (d) Fragments separated only by columns that were ALREADY gap in both the
#     parent and the child are one event -- real CNGB3 Node92->Dasypterus_ega
#     2260-2289 + 2293-2304 around an inherited 2290-2292 gap. The event length
#     is the material actually deleted (42 bp), not the 45-column span.
aln_len = 3000
gap_at = {"P": bytearray(aln_len), "A": bytearray(aln_len)}
for col in range(2290, 2293):
    gap_at["P"][col - 1] = 1
    gap_at["A"][col - 1] = 1
rows = ae.merge_contiguous_same_type_indel_events(
    [event_row(2260, 2289), event_row(2293, 2304)], gap_at=gap_at, aln_len=aln_len)
check(len(rows) == 1, f"fragments around an inherited gap merge into one event (got {len(rows)})")
if len(rows) == 1:
    check(int(rows[0]["event_length"]) == 42,
          f"the bridged event's length is the material deleted (42 bp), not the 45-column span "
          f"(got {rows[0]['event_length']})")
    check("bridging_3_inherited_gap_column" in rows[0]["breakpoint_relationships"],
          "the bridged merge records how many inherited gap columns it crossed")

# (e) The same fragments are NOT merged when the intervening columns carry real
#     retained sequence on this branch.
rows = ae.merge_contiguous_same_type_indel_events(
    [event_row(2260, 2289), event_row(2293, 2304)],
    gap_at={"P": bytearray(aln_len), "A": bytearray(aln_len)}, aln_len=aln_len)
check(len(rows) == 2,
      f"fragments separated by real retained sequence are NOT merged (got {len(rows)})")

# (f) Frame arithmetic: 1 bp + 1 bp around an inherited gap is one 2 bp
#     frameshift, not two separate frameshifts.
gap_at = {"P": bytearray(aln_len), "A": bytearray(aln_len)}
gap_at["P"][500] = gap_at["A"][500] = 1     # column 501
rows = ae.merge_contiguous_same_type_indel_events(
    [event_row(500, 500), event_row(502, 502)], gap_at=gap_at, aln_len=aln_len)
check(len(rows) == 1 and int(rows[0]["event_length"]) == 2
      and rows[0]["frame_effect"] == "frameshift",
      f"a bridged 1 bp + 1 bp pair is one 2 bp frameshift "
      f"(got {len(rows)} event(s), length {rows[0]['event_length'] if rows else 'NA'})")

if FAILURES:
    print(f"\n{len(FAILURES)} failure(s):")
    for f in FAILURES:
        print(f"  - {f}")
    raise SystemExit(1)
print("\nFunctional-shared indel / event-merging tests passed.")
