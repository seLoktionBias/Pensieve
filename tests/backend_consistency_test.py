#!/usr/bin/env python3
"""Pensieve v3.31 smoke/logic tests that do not require real external programs.

These tests are intentionally targeted at the failure modes discussed during
v3.25 review: canonical alignment preparation, breakpoint decomposition,
parsimony ties, STOP retention, runner orchestration and lesion-aware ORF
history. They do not claim biological validation of MACSE/PAML.
"""
from __future__ import annotations

import csv
import importlib.util
import os
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


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta(path):
    out = {}
    name = None
    chunks = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        out[name] = "".join(chunks)
    return out


def fasta_order(path):
    return [line[1:].split()[0] for line in Path(path).read_text().splitlines() if line.startswith(">") ]


def write_table(path, header, rows):
    with open(path, "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        w.writeheader(); w.writerows(rows)


def test_static_and_cli():
    print("static architecture and CLI")
    check((ROOT / "VERSION").read_text().strip() == "4.8", "VERSION is 4.8")
    ctl = (ROOT / "templates" / "dummy_codon_asr.ctl").read_text()
    for token, msg in [
        ("clock = 0", "codeml clock=0"),
        ("fix_blength = 0", "codeml estimates its own branch lengths"),
        ("cleandata = 0", "codeml retains ambiguous sites"),
        ("RateAncestor = 1", "codeml marginal ancestral reconstruction enabled"),
    ]:
        check(token in ctl, msg)

    runner = (ROOT / "scripts" / "run_one_gene_00_to_04.sh").read_text()
    prep_pos = runner.find('02_prepare_asr_inputs.py')
    events_pos = runner.find('03_alignment_events.py')
    backend_pos = runner.find('02_run_asr_backends.sh')
    integrate_pos = runner.find('03_integrate_asr_evidence.py')
    check(-1 not in {prep_pos, events_pos, backend_pos, integrate_pos}, "all core stages are invoked by runner")
    check(prep_pos < events_pos < backend_pos < integrate_pos,
          "runner order is alignment -> events -> ASR -> integration")
    check('candidate_indel_frameshift_events.tsv' not in runner,
          "runner no longer depends on legacy candidate-event files")
    check('paml_indelmap_asr_combined' not in runner,
          "runner no longer uses IndelMaP-projected authoritative ancestors")
    integrator = (ROOT / "scripts" / "03_integrate_asr_evidence.py").read_text()
    for obsolete in ["build_orf_table", "project_indelmap_gaps_onto_paml", "load_user_internode_registry",
                     "build_tripartition_registry", "validate_unrooted_and_rooted_tree_audit"]:
        check(f"def {obsolete}" not in integrator, f"dead legacy integration helper removed: {obsolete}")
    check("Nodes X to Y are ancestral" in integrator and "tip_count - 2" not in integrator,
          "PAML ASR completeness is driven by rst declaration, not n-2 guesswork")
    for opt in ['--min-functional-witnesses', '--tie-break', '--breakpoint-tolerance']:
        check(opt in (ROOT / "bin" / "pensieve").read_text(), f"top-level CLI forwards {opt}")
    check('--indelmap' not in (ROOT / "bin" / "pensieve").read_text(),
          "IndelMaP has been removed from the top-level CLI")

    for args in [["-h"], ["-help"], ["--help"], ["--help", "-long"]]:
        r = subprocess.run([sys.executable, str(ROOT / "bin" / "pensieve"), *args], capture_output=True, text=True)
        check(r.returncode == 0 and "Pensieve" in r.stdout, f"help form {' '.join(args)} works")
    long_help = subprocess.run([sys.executable, str(ROOT / "bin" / "pensieve"), "--help", "-long"],
                               capture_output=True, text=True).stdout
    check("Pensieve v4.8 - full manual" in long_help, "long help reports v4.8")
    check("diagnostics < alignment < events < asr < integrate < plot" in long_help,
          "long help documents actual stage order")


def test_breakpoint_decomposition():
    print("breakpoint decomposition")
    mod = load("03_alignment_events")
    intervals = [(646, 687), (646, 696), (661, 669)]
    check(mod.decompose_run((646, 696), intervals, 0) == [(646, 687), (688, 696)],
          "GUCA1B-style long run decomposes into shared core plus extension")
    check(mod.decompose_run((646, 687), intervals, 0) == [(646, 687)],
          "independent interior run does not split the shared Miniopterus core")
    check(mod.gap_runs("AA---AA-A") == [(3, 5), (8, 8)], "maximal gap runs are detected")
    # v4.4: the shared core is decomposed RECURSIVELY, so a run is broken into
    # the same atomic segments regardless of which larger run it came from. Real
    # GRK7 shape: an all-gap column 61 makes ferrumequinum a 1 bp (61,61) run,
    # most species a (61,66) run, and Molossus a longer (61,71) run. All three
    # must agree on the shared (61,61) and (62,66) atoms; Molossus then owns only
    # its genuine (67,71) extension, not a spurious undecomposed (61,66) core
    # that clusters as a Molossus-only 61-66 indel over an all-gap span.
    grk = [(61, 61), (61, 66), (61, 71)]
    check(mod.decompose_run((61, 66), grk, 0) == [(61, 61), (62, 66)],
          "the 61-66 run splits at the all-gap column into (61,61)+(62,66)")
    check(mod.decompose_run((61, 71), grk, 0) == [(61, 61), (62, 66), (67, 71)],
          "the longer 61-71 run decomposes its core consistently, not as an undecomposed (61,66)")


def test_long_deletion_does_not_manufacture_a_spurious_shared_core_event():
    # End-to-end guard for the GRK7/Molossus_molossus bug: an all-gap column
    # (col 7 here) sits at the left edge of an insertion carried by ONE tip (C,
    # residues at 8-12) while a DIFFERENT tip (D) carries a longer deletion
    # (7-15) covering the same left edge. The shared 7-12 span must never be
    # emitted as a D-only indel event -- only C's real insertion (8-12) and D's
    # own extension (13-15) are genuine events there.
    print("a long deletion sharing an all-gap left edge does not manufacture a spurious 7-12 event")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        seqs = {
            "A": "ATGAAA------AAGACC",   # gap 7-12  (majority)
            "B": "ATGAAA------AAGACC",   # gap 7-12  (majority)
            "E": "ATGAAA------AAGACC",   # gap 7-12  (majority)
            "C": "ATGAAA-CCGGGAAGACC",   # gap 7 only; residues 8-12 = its own insertion
            "D": "ATGAAA---------ACC",   # gap 7-15  (a longer deletion)
        }
        out = run_event_case(tmp, "((A:1,B:1):1,(E:1,(C:1,D:1):1):1);", seqs)
        events = tsv(out / "03_T.alignment_events.tsv")
        indels = [e for e in events if e["character_class"] == "indel"]
        span_7_12 = [e for e in indels if e["alignment_start"] == "7" and e["alignment_end"] == "12"]
        check(span_7_12 == [],
              f"no indel event spans exactly the shared all-gap-edged 7-12 region ({span_7_12})")
        d_over_insertion = [e for e in indels if e["affected_tips"] == "D"
                            and int(e["alignment_start"]) <= 12 and int(e["alignment_end"]) >= 8
                            and int(e["alignment_start"]) <= 8]
        check(d_over_insertion == [],
              f"the long-deletion tip D is not reported as its own event over C's insertion region ({d_over_insertion})")
        d_extension = [e for e in indels if e["affected_tips"] == "D"
                       and e["alignment_start"] == "13" and e["alignment_end"] == "15"]
        check(len(d_extension) == 1,
              f"D still keeps its own genuine 13-15 extension deletion ({[e['alignment_start']+'-'+e['alignment_end'] for e in indels if e['affected_tips']=='D']})")


def run_event_case(tmp: Path, tree_text: str, seqs: dict[str, str], masked_rows=None, gene="T",
                   orf_status_rows=None, extra_args=None):
    tree = tmp / "tree.nwk"; tree.write_text(tree_text + ("\n" if not tree_text.endswith("\n") else ""))
    aln = tmp / "aln.fa"; aln.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    out = tmp / "events"
    cmd = [sys.executable, str(ROOT / "scripts" / "03_alignment_events.py"),
           "--gene", gene, "--alignment", str(aln), "--tree", str(tree),
           "--outdir", str(out), "--dated", "yes", "--tie-break", "none"]
    if masked_rows is not None:
        m = tmp / "stops.tsv"
        header = ["gene", "species", "primary_alignment_start", "primary_alignment_end", "stop_codon",
                  "pseudogenizing_event_candidate", "independent_stop_candidate", "frame_shifted_at_stop"]
        write_table(m, header, masked_rows)
        cmd += ["--masked-stops", str(m)]
    if orf_status_rows is not None:
        o = tmp / "orf_status.tsv"
        write_table(o, ["gene", "species", "complete_orf"], orf_status_rows)
        cmd += ["--orf-status", str(o)]
    if extra_args:
        cmd += extra_args
    r = subprocess.run(cmd, capture_output=True, text=True)
    check(r.returncode == 0, f"event engine runs ({r.stderr[:200]})")
    return out


def test_functional_consensus_fixes_shared_indel_as_ancestral():
    # Real GUCY2F/column-2772 case: a 1 bp deletion carried by two independent
    # FUNCTIONAL (complete-ORF) lineages (A in clade (A,B), D in clade (D,E)) must
    # be reconstructed as ancestral (present at their common ancestor and lost in
    # the residue-carriers), NOT as an implausible independent deletion gained
    # convergently on each functional branch.
    print("an indel shared by independent functional lineages is fixed as ancestral, not convergent gains")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # col 5 is a gap in A and D (functional carriers); a residue in B, C, E.
        seqs = {"A": "ATGA-ACCC", "B": "ATGAAACCC", "C": "ATGAAACCC",
                "D": "ATGA-ACCC", "E": "ATGAAACCC"}
        orf = [{"gene": "T", "species": s, "complete_orf": (s in ("A", "D"))} for s in "ABCDE"]
        # WITH the rule (default min 2): the col-5 indel is fixed ancestral.
        out = run_event_case(tmp, "((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", seqs, orf_status_rows=orf,
                             extra_args=["--min-functional-witnesses", "2"])
        chars = {c["character_id"]: c for c in tsv(out / "03_T.alignment_characters.tsv")}
        col5 = [c for c in chars.values() if c["alignment_start"] == "5" and c["alignment_end"] == "5"]
        check(len(col5) == 1 and col5[0]["functional_ancestral_indel"] == "True",
              "the col-5 indel shared by two independent functional lineages is flagged functional_ancestral_indel")
        check(col5[0]["root_state"] == "present",
              "the shared functional indel is reconstructed present at the root (ancestral)")
        events = tsv(out / "03_T.alignment_events.tsv")
        cid = col5[0]["character_id"]
        gains_on_functional = [e for e in events if e["character_class"] == "indel"
                               and cid in e["event_id"] and e["event_type"] == "indel_gain"
                               and e["affected_tips"] in ("A", "D")]
        check(gains_on_functional == [],
              f"neither functional lineage carries an independent deletion of the ancestral indel ({gains_on_functional})")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        seqs = {"A": "ATGA-ACCC", "B": "ATGAAACCC", "C": "ATGAAACCC",
                "D": "ATGA-ACCC", "E": "ATGAAACCC"}
        orf = [{"gene": "T", "species": s, "complete_orf": (s in ("A", "D"))} for s in "ABCDE"]
        # Rule disabled (0): the col-5 indel is not forced ancestral.
        out = run_event_case(tmp, "((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", seqs, orf_status_rows=orf,
                             extra_args=["--min-functional-witnesses", "0"])
        chars = {c["character_id"]: c for c in tsv(out / "03_T.alignment_characters.tsv")}
        col5 = [c for c in chars.values() if c["alignment_start"] == "5" and c["alignment_end"] == "5"]
        check(len(col5) == 1 and col5[0]["functional_ancestral_indel"] == "False",
              "--min-functional-witnesses 0 disables the rule (indel not forced ancestral)")



def test_guca1b_breakpoint_history_end_to_end():
    print("GUCA1B shared-core/extension/interior history")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); length = 720; base = list("A" * length)
        def seq(gaps):
            x = base.copy()
            for start, end in gaps:
                for i in range(start - 1, end):
                    x[i] = "-"
            return "".join(x)
        seqs = {
            "Miniopterus_australis": seq([(646, 696)]),
            "Miniopterus_natalensis": seq([(646, 687)]),
            "Miniopterus_schreibersii": seq([(646, 687)]),
            "Nycteris_thebaica": seq([(661, 669)]),
            "Out1": seq([]),
            "Out2": seq([]),
        }
        tree = ("((((Miniopterus_australis:1,Miniopterus_natalensis:1):0.2,"
                "Miniopterus_schreibersii:1.2):1,Nycteris_thebaica:2.2):1,"
                "(Out1:1,Out2:1):2.2);")
        out = run_event_case(tmp, tree, seqs, gene="GUCA1B")
        rows = tsv(out / "03_GUCA1B.alignment_events.tsv")
        by_span = {(int(r["alignment_start"]), int(r["alignment_end"])): r for r in rows}
        shared = by_span.get((646, 687))
        ext = by_span.get((688, 696))
        nyc = by_span.get((661, 669))
        check(shared is not None and shared["biological_interpretation"] == "deletion"
              and shared["shared_event"] == "True" and shared["ambiguous_origin"] == "False",
              "646-687 is one confident shared Miniopterus deletion")
        check(shared is not None and set(shared["affected_tips"].split(",")) == {
              "Miniopterus_australis", "Miniopterus_natalensis", "Miniopterus_schreibersii"},
              "shared core affects exactly the three Miniopterus tips")
        check(ext is not None and ext["biological_interpretation"] == "deletion"
              and ext["origin_is_tip"] == "True" and ext["affected_tips"] == "Miniopterus_australis",
              "688-696 remains an australis-specific extension")
        check(nyc is not None and nyc["origin_is_tip"] == "True" and nyc["affected_tips"] == "Nycteris_thebaica",
              "661-669 remains an independent Nycteris event rather than splitting the Miniopterus core")


def test_contiguous_same_type_indel_fragments_are_merged_into_one_event():
    # Real bug, found by inspecting real GUCA1C events directly: Nyctalus_aviator
    # has one real, contiguous 6bp deletion (columns 400-405), but it was reported
    # as TWO separate events -- a 1bp deletion at 400, then a 5bp deletion at
    # 401-405 -- because Myotis_velifer, an unrelated tip elsewhere in the tree,
    # independently shares only the first column of that span. The breakpoint
    # decomposition correctly splits that into two characters (400 is a
    # comparable question across the whole tree, 401-405 is not), but reporting
    # them as two separate events on Nyctalus_aviator's own branch fragments one
    # real indel into an arbitrary-looking pair -- and can even manufacture fake
    # frameshift signal (1bp and 5bp are each individually not a multiple of 3,
    # even though their true combined 6bp span is in-frame).
    print("a real, contiguous, single-type indel split across characters by an unrelated tip is reported as one event")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        base = list("A" * 30)
        def seq(gaps):
            x = base.copy()
            for start, end in gaps:
                for i in range(start - 1, end):
                    x[i] = "-"
            return "".join(x)
        # Target has the real, contiguous 6bp deletion (10-15). Other, an
        # unrelated tip in a completely different clade, independently shares
        # only column 10 -- exactly the real Nyctalus_aviator/Myotis_velifer shape.
        seqs = {
            "Target": seq([(10, 15)]),
            "Sibling": seq([]),
            "Other": seq([(10, 10)]),
            "Out2": seq([]),
        }
        out = run_event_case(tmp, "((Target:1,Sibling:1):1,(Other:1,Out2:1):1);", seqs)
        rows = tsv(out / "03_T.alignment_events.tsv")
        target_rows = [r for r in rows if r["affected_tips"] == "Target" and r["shared_event"] == "False"]
        check(len(target_rows) == 1,
              f"Target's split-by-an-unrelated-tip deletion is reported as exactly one event, not several ({target_rows})")
        if len(target_rows) == 1:
            ev = target_rows[0]
            check(ev["alignment_start"] == "10" and ev["alignment_end"] == "15" and ev["event_length"] == "6",
                  f"merged event spans the real, full 10-15 deletion ({ev['alignment_start']}-{ev['alignment_end']}, {ev['event_length']}bp)")
            check(ev["biological_interpretation"] == "deletion", "merged event keeps the deletion interpretation")
            check(ev["length_mod_3"] == "0" and ev["frame_effect"] == "in_frame",
                  "the true combined 6bp span is correctly in-frame, not two spurious frameshift fragments")
            check("merged_contiguous_same_type" in ev["breakpoint_relationships"],
                  "the merge is auditable from the event row itself")
        other_rows = [r for r in rows if r["affected_tips"] == "Other" and r["shared_event"] == "False"]
        check(len(other_rows) == 1 and other_rows[0]["alignment_start"] == "10" and other_rows[0]["alignment_end"] == "10",
              f"Other's own real, independent 1bp deletion is untouched by Target's merge ({other_rows})")


def test_adjacent_different_type_indel_events_are_not_merged():
    # The flip side of the merge above, explicitly requested: an insertion
    # immediately next to a deletion on the same branch is two real, different
    # events and must stay two events, even though they are perfectly adjacent.
    print("an insertion immediately adjacent to a deletion on the same branch is never merged")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # Target loses columns 4-6 relative to everyone else (deletion) and,
        # immediately next to it, independently carries a real residue at
        # column 7 that everyone else lacks (an insertion) -- adjacent, but
        # opposite types.
        seqs = {
            "Target":  "AAA---AAAA",
            "Sibling": "AAAAAA-AAA",
            "Out1":    "AAAAAA-AAA",
            "Out2":    "AAAAAA-AAA",
        }
        out = run_event_case(tmp, "((Target:1,Sibling:1):1,(Out1:1,Out2:1):1);", seqs)
        rows = tsv(out / "03_T.alignment_events.tsv")
        target_rows = [r for r in rows if r["affected_tips"] == "Target" and r["shared_event"] == "False"]
        interps = sorted(r["biological_interpretation"] for r in target_rows)
        check(interps == ["deletion", "insertion_or_restoration"],
              f"adjacent deletion and insertion on the same branch both survive as separate events ({interps})")


def test_contiguous_shared_events_with_identical_affected_tips_are_merged():
    # Real bug, found by inspecting real CNGA3 events directly: on the
    # Node95->Node94 branch (ancestor of six Myotis/Eptesicus tips), an 18bp
    # deletion at columns 466-483 and a 3bp deletion at 484-486 -- both
    # SHARED events (shared_event=True) affecting the exact same six tips --
    # were left as two separate events even after the tip-specific merge
    # fix, because that fix explicitly excluded every shared_event=True row.
    # The same breakpoint fragmentation that splits a single tip's own event
    # also splits a shared, ancestral one; excluding shared events entirely
    # was too broad a guard. A, B and C below share a real, contiguous 10-15
    # deletion; H's own real, unrelated 10-12 deletion forces the breakpoint
    # decomposition to fragment A/B/C's run into 10-12 and 13-15 -- but both
    # fragments, reconstructed at the (A,B,C) ancestor specifically, affect
    # exactly {A,B,C}, and must be reported merged.
    print("a shared ancestral event fragmented by an unrelated tip's own partial overlap is reported as one event")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        base = list("A" * 30)
        def seq(gaps):
            x = base.copy()
            for start, end in gaps:
                for i in range(start - 1, end):
                    x[i] = "-"
            return "".join(x)
        seqs = {
            "A": seq([(10, 15)]), "B": seq([(10, 15)]), "C": seq([(10, 15)]),
            "H": seq([(10, 12)]),
            "G1": seq([]), "G2": seq([]), "G3": seq([]), "G4": seq([]),
        }
        tree = "(((A:1,B:1):0.2,C:1):1,((H:1,G1:1):1,(G2:1,(G3:1,G4:1):1):1):1);"
        out = run_event_case(tmp, tree, seqs)
        rows = tsv(out / "03_T.alignment_events.tsv")
        abc_rows = [r for r in rows if r["affected_tips"] == "A,B,C"]
        check(len(abc_rows) == 1,
              f"the (A,B,C) ancestor's split-by-H deletion is reported as exactly one shared event ({abc_rows})")
        if len(abc_rows) == 1:
            ev = abc_rows[0]
            check(ev["alignment_start"] == "10" and ev["alignment_end"] == "15" and ev["event_length"] == "6",
                  f"merged shared event spans the real, full 10-15 deletion ({ev['alignment_start']}-{ev['alignment_end']}, {ev['event_length']}bp)")
            check(ev["shared_event"] == "True", "the merged event is still correctly reported as shared")
        h_rows = [r for r in rows if r["affected_tips"] == "H"]
        check(len(h_rows) == 1 and h_rows[0]["alignment_start"] == "10" and h_rows[0]["alignment_end"] == "12",
              f"H's own real, independent 3bp deletion is untouched by the (A,B,C) merge ({h_rows})")


def test_events_with_different_affected_tips_are_not_merged_even_on_the_same_branch():
    # The other real failure mode found on the same CNGA3 branch,
    # Node95->Node94: a 45bp event affecting six tips and a 66bp event
    # affecting a completely different, non-overlapping set of 20 tips both
    # reconstruct their origin on that exact same branch (Sankoff parsimony
    # can genuinely place two unrelated characters' cheapest transition point
    # on the same branch even though their own affected-tip sets share no
    # members) -- and, separately, a contiguous pair on that same branch with
    # different affected tips is exactly the shape the affected-tips guard
    # exists for. Unit-tests the merge function directly (constructing this
    # exact origin coincidence through Sankoff reconstruction is impractical)
    # to confirm sharing only a branch, without matching affected tips, is
    # never sufficient to merge two events.
    print("two events on the same branch with different affected tips are never merged, even when contiguous")
    mod = load("03_alignment_events")
    def row(start, end, tips, event_id):
        length = end - start + 1
        return {
            "gene": "T", "event_id": event_id, "character_class": "indel",
            "event_type": "indel_gain", "biological_interpretation": "deletion",
            "alignment_start": start, "alignment_end": end, "event_length": length,
            "length_mod_3": length % 3, "frame_effect": "in_frame" if length % 3 == 0 else "frameshift",
            "origin_node": "NodeX", "origin_is_tip": False, "parent_node": "NodeY",
            "branch": "NodeY->NodeX", "shared_event": True, "n_affected_tips": len(tips.split(",")),
            "affected_tips": tips, "reversal_below_origin": False, "secondary_changes_below_origin": "NA",
            "root_state": "absent", "parsimony_score": "1", "delta_parsimony_support": "inf",
            "ambiguous_origin": False, "direction_confident": True,
            "parent_age": "2.0", "child_age": "1.0", "age_interval": "2.0-1.0",
            "n_observed_present": len(tips.split(",")), "n_observed_absent": 1, "n_unknown": 0,
            "observed_present_tips": tips, "breakpoint_relationships": "NA",
            "coordinate_system": "primary_codon_alignment",
        }
    rows = [
        row(10, 12, "A,B", "T|IND0001|gain0001"),
        row(13, 15, "C,D", "T|IND0002|gain0002"),  # contiguous, same branch, DISJOINT tips
    ]
    merged = mod.merge_contiguous_same_type_indel_events(rows)
    check(len(merged) == 2, f"two events with disjoint affected tips stay two events even when contiguous on the same branch ({merged})")


def test_overlapping_affected_tips_contiguous_events_are_clubbed():
    # v4.6, requested directly (real CNGA3 Node102->Node101): a 2128-2131 4bp
    # deletion (carried by Micronycteris_megalotis, Myotis_myotis) and a
    # contiguous 2132-2133 2bp deletion (carried by Carollia, Glyphonycteris,
    # Micronycteris_megalotis, Myotis_myotis, Trinycteris) are one 6bp in-frame
    # deletion at that branch, later partially reverted in some descendants --
    # NOT two separate frameshift deletions. Because the two fragments' affected
    # tips OVERLAP (they share carriers), that overlap is positive evidence of a
    # single common origin, so they are clubbed by position; the merged event
    # spans the full 6bp, is in-frame, and takes the union of affected tips.
    print("contiguous same-branch same-type fragments with OVERLAPPING affected tips are clubbed by position")
    mod = load("03_alignment_events")
    def row(start, end, tips, event_id):
        length = end - start + 1
        return {
            "gene": "T", "event_id": event_id, "character_class": "indel",
            "event_type": "indel_gain", "biological_interpretation": "deletion",
            "alignment_start": start, "alignment_end": end, "event_length": length,
            "length_mod_3": length % 3, "frame_effect": "in_frame" if length % 3 == 0 else "frameshift",
            "origin_node": "Node101", "origin_is_tip": False, "parent_node": "Node102",
            "branch": "Node102->Node101", "shared_event": True, "n_affected_tips": len(tips.split(",")),
            "affected_tips": tips, "reversal_below_origin": False, "secondary_changes_below_origin": "NA",
            "root_state": "absent", "parsimony_score": "1", "delta_parsimony_support": "inf",
            "ambiguous_origin": False, "direction_confident": True,
            "parent_age": "2.0", "child_age": "1.0", "age_interval": "2.0-1.0",
            "n_observed_present": len(tips.split(",")), "n_observed_absent": 1, "n_unknown": 0,
            "observed_present_tips": tips, "breakpoint_relationships": "NA",
            "coordinate_system": "primary_codon_alignment",
        }
    rows = [
        row(2128, 2131, "Micronycteris,Myotis_myotis", "T|IND0047|gain0001"),
        row(2132, 2133, "Carollia,Glyphonycteris,Micronycteris,Myotis_myotis,Trinycteris", "T|IND0048|gain0002"),
    ]
    merged = mod.merge_contiguous_same_type_indel_events(rows)
    check(len(merged) == 1, f"the overlapping contiguous 4bp+2bp pair is clubbed into one event ({len(merged)})")
    m = merged[0]
    check(int(m["alignment_start"]) == 2128 and int(m["alignment_end"]) == 2133 and int(m["event_length"]) == 6,
          f"the clubbed event spans the full 2128-2133 6bp region ({m['alignment_start']}-{m['alignment_end']})")
    check(int(m["length_mod_3"]) == 0 and m["frame_effect"] == "in_frame",
          "the clubbed 4bp+2bp deletion is correctly reported as a single 6bp IN-FRAME deletion")
    check(set(m["affected_tips"].split(",")) ==
          {"Micronycteris", "Myotis_myotis", "Carollia", "Glyphonycteris", "Trinycteris"},
          f"the clubbed event takes the union of affected tips ({m['affected_tips']})")


def test_nested_interior_event_does_not_block_merging_the_events_around_it():
    # Real bug, found by inspecting real CNGA3 events directly: Myotis_auriculus's
    # own 466-483 and 484-486 deletions were genuinely contiguous, same branch,
    # same type, same single affected tip -- and still did not merge, because a
    # THIRD event on that same branch (472-474, strictly nested inside 466-483,
    # the same "independent interior event" shape already covered by the
    # GUCA1B/Nycteris breakpoint-decomposition test) sorted in BETWEEN them by
    # start column. The merge scanned adjacency against that nested event
    # instead of against 484-486 and never saw the real contiguity at all.
    # Unit-tests the merge function directly with exactly this shape.
    print("a nested interior event on the same branch does not interrupt merging the events around it")
    mod = load("03_alignment_events")
    def row(start, end, event_id):
        length = end - start + 1
        return {
            "gene": "T", "event_id": event_id, "character_class": "indel",
            "event_type": "indel_loss", "biological_interpretation": "ambiguous_indel_change",
            "alignment_start": start, "alignment_end": end, "event_length": length,
            "length_mod_3": length % 3, "frame_effect": "in_frame" if length % 3 == 0 else "frameshift",
            "origin_node": "X", "origin_is_tip": True, "parent_node": "Y",
            "branch": "Y->X", "shared_event": False, "n_affected_tips": 1,
            "affected_tips": "X", "reversal_below_origin": False, "secondary_changes_below_origin": "NA",
            "root_state": "absent", "parsimony_score": "1", "delta_parsimony_support": "inf",
            "ambiguous_origin": False, "direction_confident": False,
            "parent_age": "2.0", "child_age": "1.0", "age_interval": "2.0-1.0",
            "n_observed_present": 1, "n_observed_absent": 1, "n_unknown": 0,
            "observed_present_tips": "X", "breakpoint_relationships": "NA",
            "coordinate_system": "primary_codon_alignment",
        }
    rows = [
        row(466, 483, "T|IND0031|loss0050"),   # outer piece 1
        row(472, 474, "T|IND0032|loss0051"),   # nested strictly inside piece 1 -- must stay separate
        row(484, 486, "T|IND0033|loss0056"),   # outer piece 2, contiguous with piece 1
    ]
    merged = mod.merge_contiguous_same_type_indel_events(rows)
    check(len(merged) == 2, f"the nested event and the merged 466-486 outer event give exactly two results ({merged})")
    spans = {(int(r["alignment_start"]), int(r["alignment_end"])) for r in merged}
    check(spans == {(466, 486), (472, 474)},
          f"the two outer pieces merge into 466-486 while the nested 472-474 event stays untouched ({spans})")


def test_parsimony_tie_and_direction():
    print("parsimony ties and insertion/deletion direction")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # 2-of-4 monophyletic gap state has two equally parsimonious polarities:
        # deletion on AB OR ancestral gap plus insertion on CD. It must stay ambiguous.
        seqs = {"A":"ATG---CCCGGG", "B":"ATG---CCCGGG",
                "C":"ATGAAACCCGGG", "D":"ATGAAACCCGGG"}
        out = run_event_case(tmp, "((A:1,B:1):1,(C:1,D:1):1);", seqs)
        rows = tsv(out / "03_T.alignment_events.tsv")
        check(rows, "tie case emits representative event rows")
        check(all(r["direction_confident"] == "False" for r in rows),
              "exact root-polarity tie never becomes direction-confident")
        check(all(r["ambiguous_origin"] == "True" for r in rows),
              "exact root-polarity tie is explicitly flagged")
        check(all(r["biological_interpretation"] == "ambiguous_indel_change" for r in rows),
              "tie is not mislabeled as a definitive deletion")

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # Three non-monophyletic gap carriers vs one residue tip strongly favors
        # ancestral GAP and a residue gain on D: an insertion/restoration, proving
        # native gap coding does not equate MACSE-style missing placeholders with deletion.
        seqs = {"A":"ATG--ACCCGGG", "B":"ATG--ACCCGGG",
                "C":"ATG--ACCCGGG", "D":"ATG-AA CCCGGG".replace(" ", "")}
        # At column 5: A/B/C have '-', D has 'A'; column 4 is gap in all tips.
        out = run_event_case(tmp, "((A:1,C:1):1,(B:1,D:1):1);", seqs)
        rows = tsv(out / "03_T.alignment_events.tsv")
        ins = [r for r in rows if r["biological_interpretation"] == "insertion_or_restoration"
               and r["direction_confident"] == "True"]
        check(bool(ins), "gap/residue context can reconstruct a confident insertion/residue gain")


def test_stop_alleles_are_separate():
    print("premature STOP event identity")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        seqs = {"A":"ATGAAATGACCC", "B":"ATGAAATAACCC",
                "C":"ATGAAACCCCCC", "D":"ATGAAACCCCCC"}
        masked = [
            {"gene":"T","species":"A","primary_alignment_start":7,"primary_alignment_end":9,
             "stop_codon":"TGA","pseudogenizing_event_candidate":True,"independent_stop_candidate":True},
            {"gene":"T","species":"B","primary_alignment_start":7,"primary_alignment_end":9,
             "stop_codon":"TAA","pseudogenizing_event_candidate":True,"independent_stop_candidate":True},
        ]
        out = run_event_case(tmp, "((A:1,B:1):1,(C:1,D:1):1);", seqs, masked)
        chars = tsv(out / "03_T.alignment_characters.tsv")
        stops = sorted(r["stop_codon"] for r in chars if r["character_class"] == "stop_mask")
        check(stops == ["TAA", "TGA"], "TAA and TGA at the same codon are separate event characters")


def test_shared_stop_not_missed_when_only_one_tip_is_registered():
    # Real bug found on a genuine HPC run (PDE6H): two sister tips (A, B) carry
    # the exact same premature-stop allele at the exact same canonical
    # alignment columns. Only A's occurrence was classified by the upstream
    # raw-sequence scan as an "independent_stop_candidate" (pseudogenizing_
    # event_candidate=True); B's occurrence fell through to a generic
    # catch-all classification (pseudogenizing_event_candidate=False) even
    # though B's canonical-alignment segment at those exact columns is the
    # identical stop codon. Before the fix, stop_mask_characters() only ever
    # marked a tip PRESENT via that per-species classification gate, so B was
    # silently treated as ABSENT despite unambiguously carrying the stop in
    # the canonical alignment -- corrupting the shared-ancestry reconstruction
    # (the true shared internal-branch origin was reported as two independent
    # tip-only gains, with no event ever attributed to the shared ancestor).
    print("shared STOP allele is not missed when only one carrier was classified as independent")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        seqs = {"A": "ATGAAATGACCC", "B": "ATGAAATGACCC",
                "C": "ATGAAACCCCCC", "D": "ATGAAACCCCCC"}
        masked = [
            {"gene": "T", "species": "A", "primary_alignment_start": 7, "primary_alignment_end": 9,
             "stop_codon": "TGA", "pseudogenizing_event_candidate": True, "independent_stop_candidate": True},
            {"gene": "T", "species": "B", "primary_alignment_start": 7, "primary_alignment_end": 9,
             "stop_codon": "TGA", "pseudogenizing_event_candidate": False, "independent_stop_candidate": False},
        ]
        out = run_event_case(tmp, "((A:1,B:1):1,(C:1,D:1):1);", seqs, masked)
        chars = tsv(out / "03_T.alignment_characters.tsv")
        stop_chars = [r for r in chars if r["character_class"] == "stop_mask"]
        check(len(stop_chars) == 1 and stop_chars[0]["n_observed_present"] == "2",
              f"both A and B are observed present for the shared STOP character ({stop_chars})")
        states = tsv(out / "03_T.alignment_character_node_states.tsv")
        b_state = next((r for r in states if r["character_id"] == "STOP0001" and r["node_label"] == "B"), None)
        check(b_state is not None and b_state["representative_state"] == "present",
              f"the unregistered-but-identical carrier B is reconstructed as present, not absent ({b_state})")
        events = tsv(out / "03_T.alignment_events.tsv")
        stop_events = [r for r in events if r["character_class"] == "stop_mask"]
        check(len(stop_events) == 1, f"one shared STOP event, not two independent ones ({stop_events})")
        ev = stop_events[0]
        check(ev["shared_event"] == "True" and ev["origin_is_tip"] == "False"
              and set(ev["affected_tips"].split(",")) == {"A", "B"},
              f"the STOP event is correctly attributed to the shared ancestor of A and B ({ev})")


def test_compensated_frameshift_stop_classification():
    print("STOP classification after compensated MACSE frame shifts")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); r0 = tmp/"r0"; r1 = tmp/"r1"; r2 = tmp/"r2"
        r0.mkdir(); r1.mkdir(); r2.mkdir()
        raw_a = "ATGAAACCCTGAGGGCCC"  # premature TGA at raw nt 10-12
        # B is the complete ORF here, so it must survive the v4.8 pre-PAML gate
        # (validate_complete_orf_alignment): its aligned row has to degap back to
        # exactly this sequence, carry no '!' and sit on whole codon cells. It is
        # therefore three codons longer than A and carries real residues opposite
        # A's frame-restoration markers, instead of the sub-codon gaps an earlier
        # version of this fixture used -- which described an alignment Pensieve now
        # (correctly) refuses to hand to PAML.
        raw_b = "ATGCAAAGGCCCCCCGGGCCC"
        (r0/"00_T.common_species.gapless_for_macse.fasta").write_text(f">A\n{raw_a}\n>B\n{raw_b}\n")
        (r0/"00_T.common_species.fasta").write_text(f">A\n{raw_a}\n>B\n{raw_b}\n")
        (r0/"00_T.common_species.tree").write_text("(A:1,B:1);\n")
        write_table(r0/"00_T.orf_status.tsv", ["species","premature_stop_count","complete_orf"], [
            {"species":"A","premature_stop_count":1,"complete_orf":False},
            {"species":"B","premature_stop_count":0,"complete_orf":True},
        ])
        write_table(r0/"00_T.orf_failures.tsv",
                    ["gene","species","failure_type","codon_position","nt_start","nt_end","codon","details"], [
            {"gene":"T","species":"A","failure_type":"premature_in_frame_stop","codon_position":4,
             "nt_start":10,"nt_end":12,"codon":"TGA","details":"raw_unaligned_sequence"},
        ])
        # Two upstream MACSE partial-codon runs have lengths 1 and 2. Their
        # cumulative frame correction is 3 mod 3 = 0, so the STOP is back in
        # phase and must remain an independent nonsense candidate.
        macse_a = "ATG!AAA!!CCCTGAGGGCCC"
        macse_b = "ATGCAAAGGCCCCCCGGGCCC"
        (r1/"01_T.macse_NT.fasta").write_text(f">A\n{macse_a}\n>B\n{macse_b}\n")
        (r1/"01_T.macse_AA.fasta").write_text(">A\nXXXXXXX\n>B\nXXXXXXX\n")
        r = subprocess.run([sys.executable, str(ROOT/"scripts"/"01_run_macse_and_extract_events.py"),
                            "--gene","T","--step00-dir",str(r0),"--outdir",str(r1)],
                           capture_output=True,text=True)
        check(r.returncode==0, f"step01 compensated-frame diagnostic runs ({r.stderr[:160]})")
        diag = tsv(r1/"01_T.macse_premature_stop_masking.tsv")
        a = [x for x in diag if x["species"]=="A"]
        check(a and a[0]["upstream_macse_marker_count"]=="2"
              and a[0]["upstream_macse_frame_correction_mod3"]=="0"
              and a[0]["frame_shifted_at_stop"]=="False",
              "two upstream MACSE corrections summing to 0 mod 3 restore frame at the STOP")
        check(a and "independent_candidate_after_compensated_frame_restoration" in a[0]["interpretation"],
              "restored-frame STOP is not mislabeled as a frameshift consequence")

        r = subprocess.run([sys.executable, str(ROOT/"scripts"/"02_prepare_asr_inputs.py"),
                            "--gene","T","--results01-dir",str(r1),"--results00-dir",str(r0),
                            "--outdir",str(r2),"--alignment-mode","perform"],capture_output=True,text=True)
        check(r.returncode==0, f"step02 accepts compensated-frame STOP ({r.stderr[:160]})")
        reg = tsv(r2/"02_T.masked_inframe_premature_stops_after_macse_correction.tsv")
        a2 = [x for x in reg if x["species"]=="A" and x["stop_codon"]=="TGA"]
        check(a2 and a2[0]["independent_stop_candidate"]=="True"
              and a2[0]["pseudogenizing_event_candidate"]=="True",
              "genuine STOP downstream of compensated frameshifts remains in the event catalogue")


def test_build_stop_registry_rejects_non_contiguous_codon_mapping():
    # v4.2 (ChatGPT review, spec 5): a raw STOP whose three nucleotides map
    # to non-consecutive canonical alignment columns (a MACSE '!'/gap wedged
    # inside the codon) must never become a pseudogenizing_event_candidate,
    # even when the frame is nominally "not shifted" at that position.
    print("build_stop_registry rejects a non-contiguous/non-3bp STOP mapping")
    mod = load("02_prepare_asr_inputs")
    # Species X: raw nt 7-9 spell TGA, but the canonical alignment has a '!'
    # wedged between raw nt 7 and 8 (columns 7,9,10 -- not contiguous, and
    # the 7..10 span is 4 columns wide, not a clean 3bp codon window).
    canonical_source = {"X": "ATGAAAT!GACCC"}
    raw_by_species = {"X": "ATGAAATGACCC"}
    raw_stops = [{
        "species": "X", "raw_codon_position": 3, "raw_nt_start": 7, "raw_nt_end": 9,
        "stop_codon": "TGA", "frame_shifted_at_stop": False,
        "upstream_macse_marker_count": 0, "upstream_macse_frame_correction_mod3": 0,
        "stop_phase_interpretation": "raw_stop_independent_candidate",
    }]
    registry = mod.build_stop_registry("T", raw_stops, canonical_source, raw_by_species, "macse_codon_alignment")
    check(len(registry) == 1, "one registry row produced")
    row = registry[0]
    check(row["mapped_columns_contiguous"] is False, f"non-contiguous mapping correctly detected ({row})")
    check(row["pseudogenizing_event_candidate"] is False and row["independent_stop_candidate"] is False,
          f"a non-contiguous/non-3bp STOP mapping is never a pseudogenizing_event_candidate ({row})")
    check("non_contiguous_or_non_3bp" in row["reason"], f"reason explains the rejection ({row['reason']})")

    # Control: the same raw STOP with a clean, contiguous alignment (no '!'
    # inside the codon) must be accepted, confirming the rejection above is
    # specifically about contiguity, not some other regression.
    clean_source = {"X": "ATGAAATGACCC"}
    clean_registry = mod.build_stop_registry("T", raw_stops, clean_source, raw_by_species, "macse_codon_alignment")
    crow = clean_registry[0]
    check(crow["mapped_columns_contiguous"] is True and crow["pseudogenizing_event_candidate"] is True,
          f"the identical raw STOP is accepted once its mapping is contiguous ({crow})")


def test_build_stop_registry_rejects_off_frame_codon_coincidence():
    # Real bug, found running real CNGA3 data (reported directly by the
    # user, confirmed and shown to them directly): a raw STOP can map to
    # THREE CONTIGUOUS canonical alignment columns that DO spell TAA/TAG/TGA
    # -- passing both the old contiguous and codon_confirmed gates -- while
    # still being pure coincidence, because those three columns do not
    # start on a real codon boundary of the alignment's own shared codon
    # frame. This happens whenever a species has an earlier synthetic
    # insertion (from an accepted conserved-block/frame correction) that is
    # itself frame-neutral in aggregate (its own total length happens to be
    # a multiple of 3) but still shifts which alignment COLUMN each later
    # raw position falls on relative to the codon grid, because the raw
    # premature-stop scanner (00_prune_and_check_orf.py) numbers codons by
    # counting from raw position 1 with no knowledge of any correction.
    #
    # Minimal real-shape reproduction: raw "ATGAAATGACCC" has a genuine raw-
    # frame stop (TGA) at raw nt 7-9, a clean codon boundary in raw terms.
    # A single synthetic 'N' inserted earlier (between raw nt 1 and 2) --
    # itself node-for-node exactly the same shape as the real CNGA3 case --
    # shifts every later raw position's alignment column by +1, so raw nt
    # 7-9 now maps to alignment columns 8-10: still contiguous, still
    # spells "TGA" verbatim (the underlying characters never change), but
    # column 8 is NOT a codon boundary (columns 1-3/4-6/7-9/10-12 are; (8-1)
    # mod 3 = 1, not 0) -- it is the tail of the real codon at columns 7-9
    # plus the head of the real codon at columns 10-12, not a real codon.
    print("build_stop_registry rejects a contiguous, textually-matching STOP that isn't on a real codon boundary")
    mod = load("02_prepare_asr_inputs")
    canonical_source = {"X": "A" + "N" + "TGAAATGACCC"}  # 13 columns; raw nt 7-9 -> columns 8-10
    raw_by_species = {"X": "ATGAAATGACCC"}
    raw_stops = [{
        "species": "X", "raw_codon_position": 3, "raw_nt_start": 7, "raw_nt_end": 9,
        "stop_codon": "TGA", "frame_shifted_at_stop": False,
        "upstream_macse_marker_count": 0, "upstream_macse_frame_correction_mod3": 0,
        "stop_phase_interpretation": "raw_stop_independent_candidate",
    }]
    registry = mod.build_stop_registry("T", raw_stops, canonical_source, raw_by_species, "macse_codon_alignment")
    row = registry[0]
    check(row["primary_alignment_start"] == 8 and row["mapped_columns_contiguous"] is True,
          f"the STOP maps to a contiguous 3-column span starting at column 8 ({row})")
    check(row["codon_frame_aligned"] is False,
          f"column 8 is correctly identified as NOT a real codon boundary ({row})")
    check(row["pseudogenizing_event_candidate"] is False and row["independent_stop_candidate"] is False,
          f"an off-frame textual coincidence is never a pseudogenizing_event_candidate, even though it "
          f"spells a stop codon verbatim ({row})")
    check("does_not_start_on_a_real_codon_boundary" in row["reason"], f"reason explains the rejection ({row['reason']})")

    # Control: the identical raw STOP without the earlier synthetic 'N' (so
    # raw nt 7-9 maps to its own natural codon boundary, columns 7-9) must
    # still be accepted -- confirming the rejection above is specifically
    # about frame alignment, not a regression in the ordinary case.
    clean_registry = mod.build_stop_registry("T", raw_stops, {"X": "ATGAAATGACCC"}, raw_by_species, "macse_codon_alignment")
    crow = clean_registry[0]
    check(crow["codon_frame_aligned"] is True and crow["pseudogenizing_event_candidate"] is True,
          f"the same STOP is accepted once it lands on a real codon boundary ({crow})")


def test_frame_dependent_stop_not_resurrected_by_coincidental_allele_match():
    # v4.2: the exact CNGA3/Phyllostomus discolor mechanism. Tip A's own raw
    # STOP occurrence at a given canonical span was classified frame-
    # dependent (a downstream consequence of an earlier, unrelated
    # frameshift) and correctly excluded from founding a character. Tip C,
    # unrelated to A's frameshift, happens to carry a genuinely independent,
    # validated nonsense mutation at the EXACT SAME canonical span/allele.
    # A must NOT be resurrected as a carrier of C's character merely because
    # A's own alignment segment there also spells the same three letters --
    # this is the precise bug that inflated CNGA3's reported event count.
    print("a frame-dependent STOP is not resurrected as a carrier via coincidental allele match")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        seqs = {"A": "ATGAAATGACCC", "B": "ATGAAACCCCCC",
                "C": "ATGAAATGACCC", "D": "ATGAAACCCCCC"}
        masked = [
            {"gene": "T", "species": "A", "primary_alignment_start": 7, "primary_alignment_end": 9,
             "stop_codon": "TGA", "pseudogenizing_event_candidate": False,
             "independent_stop_candidate": False, "frame_shifted_at_stop": True},
            {"gene": "T", "species": "C", "primary_alignment_start": 7, "primary_alignment_end": 9,
             "stop_codon": "TGA", "pseudogenizing_event_candidate": True,
             "independent_stop_candidate": True, "frame_shifted_at_stop": False},
        ]
        out = run_event_case(tmp, "((A:1,B:1):1,(C:1,D:1):1);", seqs, masked)
        chars = tsv(out / "03_T.alignment_characters.tsv")
        stop_chars = [r for r in chars if r["character_class"] == "stop_mask"]
        check(len(stop_chars) == 1 and stop_chars[0]["n_observed_present"] == "1",
              f"only C is observed present; A's frame-dependent occurrence is not resurrected "
              f"({stop_chars})")
        states = tsv(out / "03_T.alignment_character_node_states.tsv")
        a_state = next((r for r in states if r["character_id"] == "STOP0001" and r["node_label"] == "A"), None)
        check(a_state is not None and a_state["representative_state"] == "absent",
              f"A is reconstructed absent for this character despite its own segment spelling the "
              f"same codon ({a_state})")
        events = tsv(out / "03_T.alignment_events.tsv")
        stop_events = [r for r in events if r["character_class"] == "stop_mask"]
        check(len(stop_events) == 1 and stop_events[0]["affected_tips"] == "C",
              f"exactly one STOP event, attributed only to C, not A ({stop_events})")


def test_terminal_stop_stripped_from_every_input_sequence():
    # v4.0: Pensieve no longer reconstructs a terminal stop codon at
    # ancestral nodes at all (three successive designs across v3.35-v3.37
    # each fixed one real failure mode only to surface another). Instead the
    # terminal stop is removed from every input sequence up front, in
    # 00_prune_and_check_orf.py, before MACSE/MUSCLE or anything else ever
    # sees it -- there is then nothing left to reconstruct or scatter.
    print("a trailing in-frame stop codon is stripped from every sequence at input, gapless or user-curated")
    mod = load("00_prune_and_check_orf")
    check(mod.strip_terminal_stop("ATGAAACCCTGA") == "ATGAAACCC",
          "a real trailing stop is removed from a plain gapless sequence")
    check(mod.strip_terminal_stop("ATGAAACCCTAG") == "ATGAAACCC", "TAG is recognised too")
    check(mod.strip_terminal_stop("ATGAAACCCTAA") == "ATGAAACCC", "TAA is recognised too")
    check(mod.strip_terminal_stop("ATGAAACCC") == "ATGAAACCC",
          "a sequence with no terminal stop at all is left completely unchanged")
    check(mod.strip_terminal_stop("ATGAAACCCGGG") == "ATGAAACCCGGG",
          "a real terminal codon that merely isn't a stop is never removed")
    # A user-curated --alignment defined sequence may carry trailing gap
    # characters after its own real terminal stop; the 3 real stop-codon
    # bases must still be found and removed, not just the literal last 3
    # characters of the string.
    check(mod.strip_terminal_stop("ATGAAACCCTGA---") == "ATGAAACCC---",
          "trailing alignment gaps after a real terminal stop are preserved, only the 3 real bases are removed")
    check(len(mod.strip_terminal_stop("ATGAAACCCTGA")) == len("ATGAAACCCTGA") - 3,
          "exactly 3 characters are removed, never more")
    # v4.4: for --alignment defined the alignment width is authoritative and
    # must not change. preserve_width removes the terminal stop by GAPPING its
    # three real bases in place, wherever they fall (before trailing gaps),
    # instead of deleting them -- otherwise only some rows shorten and the
    # alignment desynchronises (real PDE6C/CNGB3/GUCY2F data).
    check(mod.strip_terminal_stop("ATGAAACCCTGA", preserve_width=True) == "ATGAAACCC---",
          "preserve_width gaps the terminal stop in place instead of deleting it")
    check(mod.strip_terminal_stop("ATGAAACCCTGA---", preserve_width=True) == "ATGAAACCC------",
          "preserve_width gaps the 3 real stop bases even when trailing gaps follow, keeping total width")
    check(len(mod.strip_terminal_stop("ATGAAACCCTGA", preserve_width=True)) == len("ATGAAACCCTGA"),
          "preserve_width never changes the sequence length (alignment columns are preserved)")
    check(mod.strip_terminal_stop("ATGAAACCCGGG", preserve_width=True) == "ATGAAACCCGGG",
          "preserve_width still leaves a non-stop terminal codon untouched")


def make_step_dirs(tmp: Path, gene="T"):
    r0 = tmp / "r0"; r1 = tmp / "r1"; r2 = tmp / "r2"
    r0.mkdir(); r1.mkdir(); r2.mkdir()
    (r0 / f"00_{gene}.common_species.tree").write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
    return r0, r1, r2


def test_canonical_alignment_prepare():
    print("canonical alignment preparation")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); r0, r1, r2 = make_step_dirs(tmp)
        user = {"A":"ATGAAATAACCC", "B":"ATGAAACCCCC", "C":"ATGAAACCCCCC", "D":"ATGAAACCCCCC"}
        # B is intentionally length 11 in the raw sequence: no premature STOP gate
        # should stop MACSE's frameshift-aware representation from being used.
        (r0 / "00_T.common_species.fasta").write_text("".join(f">{k}\n{v}\n" for k,v in user.items()))
        write_table(r0 / "00_T.orf_failures.tsv",
                    ["gene","species","failure_type","codon_position","nt_start","nt_end","codon","details"],
                    [{"gene":"T","species":"A","failure_type":"premature_in_frame_stop","codon_position":3,
                      "nt_start":7,"nt_end":9,"codon":"TAA","details":"raw_unaligned_sequence"},
                     {"gene":"T","species":"B","failure_type":"length_not_multiple_of_3","codon_position":"NA",
                      "nt_start":"NA","nt_end":"NA","codon":"NA","details":"length=11"}])
        macse = {"A":"ATGAAATAACCC", "B":"ATGAAA!CCCCC", "C":"ATGAAACCCCCC", "D":"ATGAAACCCCCC"}
        (r1 / "01_T.macse_NT.fasta").write_text("".join(f">{k}\n{v}\n" for k,v in macse.items()))
        write_table(r1 / "01_T.macse_premature_stop_masking.tsv",
                    ["species","raw_nt_start","raw_nt_end","masked_by_upstream_macse_frameshift_marker"],
                    [{"species":"A","raw_nt_start":7,"raw_nt_end":9,
                      "masked_by_upstream_macse_frameshift_marker":False}])
        r = subprocess.run([sys.executable, str(ROOT/"scripts"/"02_prepare_asr_inputs.py"),
                            "--gene","T","--results01-dir",str(r1),"--results00-dir",str(r0),
                            "--outdir",str(r2),"--alignment-mode","perform"], capture_output=True, text=True)
        check(r.returncode == 0, f"perform-mode preparation runs ({r.stderr[:120]})")
        native = fasta(r2 / "02_T.primary_codon_alignment_native.fasta")
        safe = fasta(r2 / "02_T.primary_codon_alignment.fasta")
        check(native["B"][6] == "-" and safe["B"][6] == "N",
              "frameshift-only taxon is retained through MACSE placeholder conversion without a STOP gate")
        check(native["A"][6:9] == "TAA" and safe["A"][6:9] == "NNN",
              "raw premature STOP survives in native coordinates and is masked only in PAML view")
        stops = tsv(r2 / "02_T.masked_inframe_premature_stops_after_macse_correction.tsv")
        a = [x for x in stops if x["species"] == "A" and x["stop_codon"] == "TAA"]
        check(a and a[0]["primary_alignment_start"] == "7" and a[0]["independent_stop_candidate"] == "True",
              "premature STOP registry has canonical coordinates before masking")
        check((r2 / "02_T.codon_for_paml.phy").stat().st_size > 0,
              "PAML PHYLIP input is actually created")

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); r0, r1, r2 = make_step_dirs(tmp)
        defined = {"A":"ATGAAA---CCC", "B":"ATGAAATTTCCC", "C":"ATGAAATTTCCC", "D":"ATGAAATTTCCC"}
        (r0 / "00_T.common_species.fasta").write_text("".join(f">{k}\n{v}\n" for k,v in defined.items()))
        write_table(r0 / "00_T.orf_failures.tsv",
                    ["gene","species","failure_type","codon_position","nt_start","nt_end","codon","details"], [])
        # Deliberately different MACSE alignment. Defined mode must ignore its columns.
        (r1 / "01_T.macse_NT.fasta").write_text("".join(f">{k}\nATGAAACCCCCC\n" for k in defined))
        write_table(r1 / "01_T.macse_premature_stop_masking.tsv",
                    ["species","raw_nt_start","raw_nt_end","masked_by_upstream_macse_frameshift_marker"], [])
        r = subprocess.run([sys.executable, str(ROOT/"scripts"/"02_prepare_asr_inputs.py"),
                            "--gene","T","--results01-dir",str(r1),"--results00-dir",str(r0),
                            "--outdir",str(r2),"--alignment-mode","defined"], capture_output=True, text=True)
        check(r.returncode == 0, f"defined-mode preparation runs ({r.stderr[:120]})")
        native = fasta(r2 / "02_T.primary_codon_alignment_native.fasta")
        check(native == defined, "defined-mode user alignment columns are preserved exactly per sequence")
        prov = tsv(r2 / "02_T.alignment_provenance.tsv")[0]
        check(prov["defined_alignment_columns_modified"] == "False" and prov["second_alignment_performed"] == "False",
              "defined mode records that no column insertion/re-alignment occurred")


def test_defined_alignment_direct_codon_stop_scan_without_macse():
    # v4.4: --alignment defined detects premature stops by reading the codon
    # alignment directly (no MACSE, no raw->alignment remap). This test does
    # NOT create any 01_*.macse_* file at all, proving step 02 defined mode
    # never depends on MACSE output. C carries an internal in-frame TGA at
    # codon 3 (columns 7-9) plus a natural terminal TAA at codon 4; only the
    # internal one may found an independent nonsense event, and both stops must
    # be masked to NNN for the codeml-safe view while native retains them.
    print("--alignment defined detects premature stops by direct codon-frame scan, no MACSE required")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); r0, r1, r2 = make_step_dirs(tmp)
        defined = {"A":"ATGAAAGGGCCC", "B":"ATGAAAGGGCCC", "C":"ATGAAATGATAA", "D":"ATGAAAGGGCCC"}
        (r0 / "00_T.common_species.fasta").write_text("".join(f">{k}\n{v}\n" for k,v in defined.items()))
        write_table(r0 / "00_T.orf_failures.tsv",
                    ["gene","species","failure_type","codon_position","nt_start","nt_end","codon","details"], [])
        r = subprocess.run([sys.executable, str(ROOT/"scripts"/"02_prepare_asr_inputs.py"),
                            "--gene","T","--results01-dir",str(r1),"--results00-dir",str(r0),
                            "--outdir",str(r2),"--alignment-mode","defined"], capture_output=True, text=True)
        check(r.returncode == 0, f"defined-mode scan runs with no MACSE file present ({r.stderr[:160]})")
        native = fasta(r2 / "02_T.primary_codon_alignment_native.fasta")
        safe = fasta(r2 / "02_T.primary_codon_alignment.fasta")
        check(native["C"][6:9] == "TGA", "native view retains C's internal stop allele unchanged")
        check(safe["C"][6:9] == "NNN" and safe["C"][9:12] == "NNN",
              "PAML-safe view masks BOTH C's internal stop and its terminal stop to NNN")
        # no stop codons whatsoever survive into the codeml-safe alignment
        surviving = [safe[s][i:i+3] for s in safe for i in range(0, len(safe[s]), 3)
                     if safe[s][i:i+3] in {"TAA","TAG","TGA"}]
        check(surviving == [], f"no stop codon survives in the PAML-safe alignment ({surviving})")
        stops = tsv(r2 / "02_T.masked_inframe_premature_stops_after_macse_correction.tsv")
        internal = [x for x in stops if x["species"]=="C" and x["primary_alignment_start"]=="7"]
        terminal = [x for x in stops if x["species"]=="C" and x["primary_alignment_start"]=="10"]
        check(len(internal)==1 and internal[0]["stop_codon"]=="TGA"
              and internal[0]["pseudogenizing_event_candidate"]=="True"
              and internal[0]["independent_stop_candidate"]=="True",
              "the internal TGA is registered as an independent premature-stop candidate at its exact columns")
        check(len(terminal)==1 and terminal[0]["terminal_codon"]=="True"
              and terminal[0]["pseudogenizing_event_candidate"]=="False",
              "the terminal stop is recorded but excluded from founding an event")
        prov = tsv(r2 / "02_T.alignment_provenance.tsv")[0]
        check(prov["aligner"]=="NA" and prov["premature_stop_detection"].startswith("direct_codon_frame_scan"),
              "provenance records the MACSE-free direct codon-frame stop detection")


def test_prepare_defined_alignment_frame_fix_and_masking():
    # v4.4: prepare_defined_alignment.py makes a curated alignment ready for
    # --alignment defined using ONLY mechanical, data-preserving operations:
    # mask flagged assembly gaps ('-' -> 'N') and trim leading/trailing all-gap
    # padding so the length becomes a codon multiple. Internal all-gap columns
    # (real codon-frame insertion columns) and every real base are untouched.
    print("prepare_defined_alignment: trailing all-gap trim to codon frame + assembly-gap masking")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # length 11: internal all-gap column at index 3 (kept), trailing all-gap
        # columns 9-10 (trimmed -> length 9, a codon multiple). A carries a
        # per-species gap at index 6 that a mask region rewrites to N.
        aln = {"A":"ATG-AA-CC--", "B":"ATG-AATCC--", "C":"ATG-TATCC--"}
        (tmp/"in.fa").write_text("".join(f">{k}\n{v}\n" for k,v in aln.items()))
        (tmp/"mask.txt").write_text("A\t7\t7\n")  # 1-based col 7 = A's gap at index 6
        r = subprocess.run([sys.executable, str(ROOT/"scripts"/"prepare_defined_alignment.py"),
                            "--in", str(tmp/"in.fa"), "--mask-regions", str(tmp/"mask.txt"),
                            "--out", str(tmp/"out.fa")], capture_output=True, text=True)
        check(r.returncode == 0, f"preprocessing runs ({r.stderr[:160]})")
        out = fasta(tmp/"out.fa")
        L = {len(v) for v in out.values()}
        check(L == {9}, f"output is trimmed to a single codon-multiple length ({L})")
        check(out["A"][6] == "N", "a flagged assembly gap is rewritten to N")
        check(out["B"] == "ATG-AATCC" and out["C"] == "ATG-TATCC",
              "real bases and the internal all-gap column are preserved untouched")
        check(all(v[3] == "-" for v in out.values()),
              "the internal all-gap codon-frame column is retained, not removed")
    # Refuses (does not mangle) when the residual is genuinely out of frame.
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        offframe = {"A":"ATGAACC", "B":"ATGAACC"}  # length 7, no all-gap padding to trim
        (tmp/"in.fa").write_text("".join(f">{k}\n{v}\n" for k,v in offframe.items()))
        r = subprocess.run([sys.executable, str(ROOT/"scripts"/"prepare_defined_alignment.py"),
                            "--in", str(tmp/"in.fa"), "--out", str(tmp/"out.fa")], capture_output=True, text=True)
        check(r.returncode != 0 and not (tmp/"out.fa").exists(),
              "a genuinely off-frame alignment is refused, not silently trimmed")


def test_phylogenetic_sequence_ordering():
    print("phylogenetic ordering of species and internodes")
    with tempfile.TemporaryDirectory() as td:
        tmp=Path(td); r0=tmp/"r0"; r1=tmp/"r1"; r2=tmp/"r2"; r0.mkdir(); r1.mkdir(); r2.mkdir()
        # Deliberately non-alphabetical left-to-right tree order: C,A,D,B.
        tree=tmp/"tree.nwk"; tree.write_text("((C:1,A:1):1,(D:1,B:1):1);\n")
        fa=tmp/"T.fa"; seq="ATGAAACCCGGG"
        fa.write_text("".join(f">{x}\n{seq}\n" for x in ["A","B","C","D"]))
        r=subprocess.run([sys.executable,str(ROOT/"scripts"/"00_prune_and_check_orf.py"),"--gene","T",
                          "--fasta",str(fa),"--tree",str(tree),"--outdir",str(r0)],capture_output=True,text=True)
        check(r.returncode==0,f"step00 phylogenetic ordering runs ({r.stderr[:120]})")
        # Requested directly: the exported order must match how the tree is
        # actually rendered top-to-bottom by 05_plot_events.R, not plain
        # left-to-right preorder (which reads bottom-to-top next to the plot).
        check(fasta_order(r0/"00_T.common_species.fasta")==["B","D","A","C"],
              "step00 FASTA follows the plot's top-to-bottom tip order rather than alphabet")

        # Pretend MACSE returned alphabetical order; step01 must rewrite its files
        # to the rooted-tree tip order even when reusing existing MACSE outputs.
        (r1/"01_T.macse_NT.fasta").write_text("".join(f">{x}\n{seq}\n" for x in ["A","B","C","D"]))
        (r1/"01_T.macse_AA.fasta").write_text("".join(f">{x}\nXXXX\n" for x in ["A","B","C","D"]))
        r=subprocess.run([sys.executable,str(ROOT/"scripts"/"01_run_macse_and_extract_events.py"),"--gene","T",
                          "--step00-dir",str(r0),"--outdir",str(r1)],capture_output=True,text=True)
        check(r.returncode==0,f"step01 reorders reused MACSE output ({r.stderr[:120]})")
        check(fasta_order(r1/"01_T.macse_NT.fasta")==["B","D","A","C"],
              "MACSE NT MSA is rewritten in the plot's top-to-bottom tip order")
        check(fasta_order(r1/"01_T.macse_AA.fasta")==["B","D","A","C"],
              "MACSE AA MSA is rewritten in the plot's top-to-bottom tip order")

        # Step02 final canonical and PAML-safe alignments must retain same order.
        write_table(r1/"01_T.macse_premature_stop_masking.tsv",
                    ["species","raw_nt_start","raw_nt_end","masked_by_upstream_macse_frameshift_marker"],[])
        r=subprocess.run([sys.executable,str(ROOT/"scripts"/"02_prepare_asr_inputs.py"),"--gene","T",
                          "--results01-dir",str(r1),"--results00-dir",str(r0),"--outdir",str(r2),
                          "--alignment-mode","perform"],capture_output=True,text=True)
        check(r.returncode==0,f"step02 phylogenetic ordering runs ({r.stderr[:120]})")
        check(fasta_order(r2/"02_T.primary_codon_alignment_native.fasta")==["B","D","A","C"],
              "canonical alignment is in the plot's top-to-bottom tip order")

    # Combined tip+internode order: root is unavailable under clock=0, then each
    # internal ancestor immediately precedes its descendant clade.
    with tempfile.TemporaryDirectory() as td:
        tmp=Path(td); out=tmp/"out"; out.mkdir()
        pens=tmp/"pens.nwk"; pens.write_text("((C:1,A:1)NodeCA:1,(D:1,B:1)NodeDB:1)UserRoot;\n")
        paml=tmp/"paml.nwk"; paml.write_text("((C:1,A:1)PCA:1,(D:1,B:1)PDB:1)PRoot;\n")
        asr=tmp/"asr.fa"; asr.write_text(">PCA\nATGAAACCCGGG\n>PDB\nATGAAACCCGGG\n")
        tips=tmp/"tips.fa"; tips.write_text("".join(f">{x}\nATGAAACCCGGG\n" for x in ["A","B","C","D"]))
        chars=tmp/"chars.tsv"; write_table(chars,["character_id","character_class","alignment_start","alignment_end","length_mod_3","stop_codon"],[])
        states=tmp/"states.tsv"; write_table(states,["character_id","node_label","state"],[])
        events=tmp/"events.tsv"; write_table(events,["event_id","branch","character_class","biological_interpretation","length_mod_3","direction_confident"],[])
        frame=tmp/"frame.tsv"; write_table(frame,["node_label","frame_currently_shifted","structural_state_ambiguous"],[])
        r=subprocess.run([sys.executable,str(ROOT/"scripts"/"04_ancestral_orf_walk.py"),"--gene","T",
                          "--paml-tree",str(paml),"--pensieve-tree",str(pens),"--ancestral-fasta",str(asr),
                          "--tip-alignment",str(tips),"--character-states",str(states),"--characters",str(chars),
                          "--events",str(events),"--frame-arithmetic",str(frame),"--outdir",str(out)],
                         capture_output=True,text=True)
        check(r.returncode==0,f"combined phylogenetic MSA construction runs ({r.stderr[:120]})")
        # Requested directly: ape::plot.phylo (05_plot_events.R) draws the FIRST
        # tip in normal left-to-right preorder at the BOTTOM of the figure, so a
        # file written in that order reads bottom-to-top next to the plot --
        # backwards from how a person reads a file top to bottom. The exported
        # order is therefore reverse-preorder (each node's own children visited
        # right to left, root still before its descendants), which reads
        # top-to-bottom in exactly the same order as the rendered tree.
        expected=["NodeDB","B","D","NodeCA","A","C"]
        check(fasta_order(out/"03_T.phylogenetic_msa.fa")==expected,
              "combined species+internode MSA follows the plot's top-to-bottom order")
        order_rows=tsv(out/"03_T.phylogenetic_sequence_order.tsv")
        check([x["node_label"] for x in order_rows]==["UserRoot","NodeDB","B","D","NodeCA","A","C"],
              "order audit table records root and all nodes in the plot's top-to-bottom order")
        check(order_rows[0]["sequence_available"]=="False" and order_rows[0]["sequence_fasta_rank"]=="NA",
              "clock=0 root is documented but not fabricated in the MSA")


def test_runner_orchestration_through_events():
    print("runner orchestration through events (mock MACSE, no real external programs)")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); bindir = tmp / "bin"; bindir.mkdir(); work = tmp / "work"; work.mkdir()
        mock = bindir / "macse"
        mock.write_text(r'''#!/usr/bin/env python3
import sys
from pathlib import Path
args=sys.argv[1:]
def val(k): return args[args.index(k)+1]
seq=Path(val('-seq')); outnt=Path(val('-out_NT')); outaa=Path(val('-out_AA'))
records=[]; name=None; chunks=[]
for line in seq.read_text().splitlines():
    if line.startswith('>'):
        if name is not None: records.append((name,''.join(chunks)))
        name=line[1:].split()[0]; chunks=[]
    else: chunks.append(line.strip())
if name is not None: records.append((name,''.join(chunks)))
with outnt.open('w') as n, outaa.open('w') as a:
    for name,s in records:
        # A/B get a whole-codon gap INSERTED (not real residues overwritten):
        # every sequence here is a complete ORF, so the v4.8 pre-PAML gate
        # requires each row to degap back to its exact input and to sit on whole
        # codon cells. C/D are padded to the same width.
        s = s[:6]+'---'+s[6:] if name in {'A','B'} else s+'---'
        n.write(f'>{name}\n{s}\n')
        a.write(f'>{name}\n'+('X'*(len(s)//3))+'\n')
''')
        mock.chmod(0o755)
        fa = tmp / "T.fasta"
        seq = "ATGCCCAAAGGGTTTCCC"
        fa.write_text("".join(f">{x}\n{seq}\n" for x in "ABCD"))
        tree = tmp / "tree.nwk"; tree.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
        env = os.environ.copy(); env["PATH"] = f"{bindir}:{env['PATH']}"
        r = subprocess.run([str(ROOT/"scripts"/"run_one_gene_00_to_04.sh"),
                            "--gene","T","--fasta",str(fa),"--tree",str(tree),"--workdir",str(work),
                            "--run_up_to","events","--alignment","perform"],
                           capture_output=True, text=True, env=env)
        check(r.returncode == 0, f"runner reaches events without missing step02 file ({r.stderr[-180:]})")
        check((work/"results_02/T/02_T.codon_for_paml.phy").stat().st_size > 0,
              "runner created the file that v3.25 previously failed to create")
        check((work/"results_03/T/03_T.alignment_events.tsv").stat().st_size > 0,
              "runner reaches and completes event reconstruction")



def test_runner_orchestration_through_integrate():
    print("runner orchestration through integrate (mock MACSE + mock codeml)")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); bindir = tmp / "bin"; bindir.mkdir(); work = tmp / "work"; work.mkdir()

        mock_macse = bindir / "macse"
        mock_macse.write_text(r'''#!/usr/bin/env python3
import sys
from pathlib import Path
args=sys.argv[1:]
def val(k): return args[args.index(k)+1]
seq=Path(val('-seq')); outnt=Path(val('-out_NT')); outaa=Path(val('-out_AA'))
records=[]; name=None; chunks=[]
for line in seq.read_text().splitlines():
    if line.startswith('>'):
        if name is not None: records.append((name,''.join(chunks)))
        name=line[1:].split()[0]; chunks=[]
    else: chunks.append(line.strip())
if name is not None: records.append((name,''.join(chunks)))
with outnt.open('w') as n, outaa.open('w') as a:
    for name,s in records:
        # A/B get a whole-codon gap INSERTED (not real residues overwritten):
        # every sequence here is a complete ORF, so the v4.8 pre-PAML gate
        # requires each row to degap back to its exact input and to sit on whole
        # codon cells. C/D are padded to the same width.
        s = s[:6]+'---'+s[6:] if name in {'A','B'} else s+'---'
        n.write(f'>{name}\n{s}\n')
        a.write(f'>{name}\n'+('X'*(len(s)//3))+'\n')
''')
        mock_macse.chmod(0o755)

        mock_codeml = bindir / "codeml"
        mock_codeml.write_text(r'''#!/usr/bin/env python3
from pathlib import Path
# 21 columns: the mock MACSE above INSERTS a whole-codon gap rather than
# overwriting residues (required by the v4.8 pre-PAML complete-ORF gate), so the
# canonical alignment is 21 wide and every marginal ancestral sequence must be too.
seq='ATGCCCAAAGGGTTTCCCAAA'
Path('rst').write_text(f"""tree with node labels for Rod Page's TreeView
((1_A:0.1,2_B:0.1)6:0.1,(3_C:0.1,4_D:0.1)7:0.1)5;
Nodes 5 to 7 are ancestral
Unreliable at sites with alignment gaps

(1) Marginal reconstruction of ancestral sequences
node #5 {seq}
node #6 {seq}
node #7 {seq}
Overall accuracy of the reconstruction
""")
Path('codon_asr.out').write_text('mock codeml output\n')
''')
        mock_codeml.chmod(0o755)

        fa = tmp / "T.fasta"
        seq = "ATGCCCAAAGGGTTTCCC"
        fa.write_text("".join(f">{x}\n{seq}\n" for x in "ABCD"))
        tree = tmp / "tree.nwk"; tree.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
        env = os.environ.copy(); env["PATH"] = f"{bindir}:{env['PATH']}"
        r = subprocess.run([str(ROOT/"scripts"/"run_one_gene_00_to_04.sh"),
                            "--gene","T","--fasta",str(fa),"--tree",str(tree),"--workdir",str(work),
                            "--run_up_to","integrate","--alignment","perform"],
                           capture_output=True, text=True, env=env)
        check(r.returncode == 0, f"runner reaches integration without v3.25 wiring failures ({r.stderr[-220:]})")
        expected = [
            (work/"results_02/T/02_T.codon_for_paml.phy", "runner creates PAML PHYLIP before invoking the backend"),
            (work/"results_02/T/paml_codon_asr/rst", "runner invokes codeml after ASR inputs exist"),
            (work/"results_03/T/03_T.alignment_events.tsv", "runner reconstructs events before ancestral integration"),
            (work/"results_03/T/03_T.ancestral_integrated_alignment.fa", "runner completes Pensieve-state overlay onto PAML scaffold"),
            (work/"results_03/T/04_T.ancestral_orf_walk.tsv", "runner completes lesion-aware ORF history"),
            (work/"final_results/T/important_output/T.ancestral_integrated_alignment.fa", "integrated core output is copied into final_results"),
            (work/"final_results/T/important_output/T.phylogenetic_msa.fasta", "combined phylogenetic tip+internode MSA is copied into final_results"),
            (work/"final_results/T/important_output/T.phylogenetic_sequence_order.tsv", "phylogenetic sequence-order audit is copied into final_results"),
            (work/"final_results/T/important_output/T.phylogenetic_tip_order.tsv", "phylogenetic tip-order audit is copied into final_results"),
        ]
        for path, message in expected:
            check(path.exists() and path.stat().st_size > 0, message)

def test_mock_paml_integration_and_root_policy():
    print("mock PAML integration and root policy")
    with tempfile.TemporaryDirectory() as td:
        tmp=Path(td); r2=tmp/"r2"; r3=tmp/"r3"; (r2/"paml_codon_asr").mkdir(parents=True); r3.mkdir()
        (r2/"02_T.codon_for_paml.phy").write_text("4 12\nA ATGAAACCCGGG\nB ATGAAACCCGGG\nC ATGAAACCCGGG\nD ATGAAACCCGGG\n")
        (r2/"02_T.tree_for_asr.nwk").write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
        aln=tmp/"aln.fa"; aln.write_text("".join(f">{x}\nATGAAACCCGGG\n" for x in "ABCD"))
        er=subprocess.run([sys.executable,str(ROOT/"scripts"/"03_alignment_events.py"),"--gene","T",
                           "--alignment",str(aln),"--tree",str(r2/"02_T.tree_for_asr.nwk"),"--outdir",str(r3)],
                          capture_output=True,text=True)
        check(er.returncode==0,"event tree prepared for PAML mapping test")
        (r2/"paml_codon_asr/rst").write_text('''tree with node labels for Rod Page's TreeView
((1_A:0.1,2_B:0.1)6:0.1,(3_C:0.1,4_D:0.1)7:0.1)5;
Nodes 5 to 7 are ancestral
Unreliable at sites with alignment gaps

(1) Marginal reconstruction of ancestral sequences
node #5 ATGAAACCCGGG
node #6 ATGAAACCCGGG
node #7 ATGAAACCCGGG
Overall accuracy of the reconstruction
''')
        ir=subprocess.run([sys.executable,str(ROOT/"scripts"/"03_integrate_asr_evidence.py"),"--gene","T",
                           "--results02-dir",str(r2),"--results00-dir",str(tmp/"r0"),"--outdir",str(r3),
                           "--on-missing-root-sequence","warn"],capture_output=True,text=True)
        check(ir.returncode==0,f"mock PAML marginal sequences map to Pensieve nodes ({ir.stderr[:120]})")
        cross=tsv(r3/"03_T.internode_label_crosswalk.tsv")
        root=[x for x in cross if x["is_user_tree_root"]=="True"]
        check(root and root[0]["paml_node_label"]=="NA" and root[0]["mapping_status"]=="NO_PAML_SEQUENCE_EXPECTED",
              "biological root is not assigned a fabricated clock=0 PAML sequence")
        all_declared=fasta(r3/"03_T.paml_marginal_asr_all_declared_nodes.fa")
        biological=fasta(r3/"03_T.paml_marginal_asr.fa")
        check(set(all_declared)=={"PAML_Node5","PAML_Node6","PAML_Node7"},
              "all rst-declared PAML marginal nodes are retained for audit")
        check(set(biological)=={"PAML_Node6","PAML_Node7"},
              "serialization/root-only PAML node is excluded from biological non-root ASR")
        orf=subprocess.run([sys.executable,str(ROOT/"scripts"/"04_ancestral_orf_walk.py"),"--gene","T",
                            "--paml-tree",str(r3/"03_T.paml_labeled_reporting_tree.nwk"),
                            "--pensieve-tree",str(r3/"03_T.pensieve_labelled_dated_tree.nwk"),
                            "--ancestral-fasta",str(r3/"03_T.paml_marginal_asr.fa"),"--tip-alignment",str(aln),
                            "--character-states",str(r3/"03_T.alignment_character_node_states.tsv"),
                            "--characters",str(r3/"03_T.alignment_characters.tsv"),
                            "--events",str(r3/"03_T.alignment_events.tsv"),
                            "--frame-arithmetic",str(r3/"03_T.frame_arithmetic_by_node.tsv"),"--outdir",str(r3)],
                           capture_output=True,text=True)
        check(orf.returncode==0,f"native ancestral ORF integration runs ({orf.stderr[:120]})")
        root_status=[x for x in tsv(r3/"04_T.ancestral_orf_walk.tsv") if x["node_type"]=="root"]
        check(root_status and root_status[0]["coding_status"]=="unavailable",
              "root ORF status remains unavailable rather than invented")


def test_partial_orf_status_distinguishes_truncation_from_pseudogenization():
    # v4.4: ORF status is three-way. A complete ORF is intact; an internal
    # premature stop or an internal (non-terminal) frameshift is disrupted
    # (pseudogenized); an incomplete sequence with NEITHER is "partial" --
    # incomplete is not evidence of pseudogenization (real GRK7 Cynopterus case).
    print("ORF status: incomplete-but-clean is partial, internal frameshift/stop is disrupted")
    mod = load("04_ancestral_orf_walk")
    check(mod.orf_check("ATGAAACCCTAA")["coding_status"] == "intact",
          "a complete ORF (ATG, mod3, no internal stop) is intact")
    check(mod.orf_check("ATGAAACC----")["coding_status"] == "partial",
          "a 3'-truncated sequence with no internal frameshift/stop is partial")
    check(mod.orf_check("GTGAAACCCGGG")["coding_status"] == "partial",
          "a missing start codon alone is partial, not disrupted")
    check(mod.orf_check("ATGTAACCCGGG")["coding_status"] == "disrupted",
          "an internal premature stop is disrupted (pseudogenized)")
    check(mod.orf_check("ATGAA-CCCGG")["coding_status"] == "disrupted",
          "an internal non-mod3 gap (internal frameshift) is disrupted")
    check(mod.orf_check("ATGA-AACCCGGG")["coding_status"] == "intact",
          "a 1bp gap whose gap-stripped sequence is still a clean ORF stays intact (compensated/artifact)")
    # has_internal_frameshift: only NON-terminal, non-mod3 gap runs count.
    check(mod.has_internal_frameshift("--ATGAAACCC--") is False,
          "leading/trailing gaps are truncation, never an internal frameshift")
    check(mod.has_internal_frameshift("ATG-AACCC") is True,
          "an internal 1bp gap is an internal frameshift")
    check(mod.has_internal_frameshift("ATG---AACCC") is False,
          "an internal 3bp gap removes a whole codon and preserves frame")
    # step 03 terminal-incompleteness helper
    ev = load("03_alignment_events")
    occ = {"X": (1, 100)}
    check(ev.event_is_terminal_incompleteness(occ, ["X"], 90, 100) is True,
          "a 3'-terminal indel is flagged as terminal incompleteness")
    check(ev.event_is_terminal_incompleteness(occ, ["X"], 1, 10) is True,
          "a 5'-terminal indel is flagged as terminal incompleteness")
    check(ev.event_is_terminal_incompleteness(occ, ["X"], 40, 50) is False,
          "an interior indel is not terminal incompleteness")


def test_root_adjacent_uncertainty_is_explicit():
    print("root-adjacent substitution-only disruption is explicit, not generic uncertain")
    with tempfile.TemporaryDirectory() as td:
        tmp=Path(td); out=tmp/"o"; out.mkdir()
        pens=tmp/"pens.nwk"; pens.write_text("((A:1,B:1)NodeAB:1,(C:1,D:1)NodeCD:1)UserRoot;\n")
        paml=tmp/"paml.nwk"; paml.write_text("((A:1,B:1)PAB:1,(C:1,D:1)PCD:1)PRoot;\n")
        # PAB carries an internal premature stop (a substitution-only disruption)
        # while PCD is intact. There is no explicit indel/stop event character, so
        # the UserRoot->NodeAB first-loss placement must remain unresolved. (A mere
        # missing ATG with no internal frameshift/stop is now "partial", not
        # disrupted, so a genuine internal stop is used to exercise this pathway.)
        asr=tmp/"asr.fa"; asr.write_text(">PAB\nATGTAACCCGGG\n>PCD\nATGAAACCCGGG\n")
        tips=tmp/"tips.fa"; tips.write_text(">A\nATGTAACCCGGG\n>B\nATGTAACCCGGG\n>C\nATGAAACCCGGG\n>D\nATGAAACCCGGG\n")
        chars=tmp/"chars.tsv"; write_table(chars,["character_id","character_class","alignment_start","alignment_end","length_mod_3","stop_codon"],[])
        states=tmp/"states.tsv"; write_table(states,["character_id","node_label","state"],[])
        events=tmp/"events.tsv"; write_table(events,["event_id","branch","character_class","biological_interpretation","length_mod_3","direction_confident"],[])
        frame=tmp/"frame.tsv"; write_table(frame,["node_label","frame_currently_shifted","structural_state_ambiguous"],[])
        r=subprocess.run([sys.executable,str(ROOT/"scripts"/"04_ancestral_orf_walk.py"),"--gene","T",
                          "--paml-tree",str(paml),"--pensieve-tree",str(pens),"--ancestral-fasta",str(asr),
                          "--tip-alignment",str(tips),"--character-states",str(states),"--characters",str(chars),
                          "--events",str(events),"--frame-arithmetic",str(frame),"--outdir",str(out)],
                         capture_output=True,text=True)
        check(r.returncode==0, f"root-adjacent uncertainty case runs ({r.stderr[:160]})")
        br={x["branch"]:x for x in tsv(out/"04_T.orf_transitions_by_branch.tsv")}
        row=br["UserRoot->NodeAB"]
        check(row["orf_transition"]=="root_adjacent_disruption_first_loss_unresolved",
              "root-adjacent sequence disruption gets a specific first-loss-unresolved label")
        check(row["uncertainty_reason"]=="no_distinct_clock0_biological_root_sequence",
              "root-adjacent asymmetry is explicitly explained in the output")
        check(row["known_pseudogenic_history"]=="True",
              "disrupted child establishes sticky pseudogenic history without inventing a root sequence")


def test_uncertain_root_resolves_to_intact_once_descendants_are_confidently_intact():
    # Real bug found by reviewing genuine HPC runs (GUCA1B/GUCA1C/CNGA3): the
    # biological root never gets a fabricated sequence, so its own
    # pseudogenic-history state starts "unresolved" for every gene by design.
    # Before this fix, once a node's history was unresolved it could only ever
    # become "confirmed disabling" -- there was no path back to "confidently
    # clean", so "history_uncertain_no_confirmed_event" (plotted as an
    # unbroken chain of amber/orange branches) silently propagated down
    # EVERY descendant branch forever, even when a node's own reconstructed
    # sequence was a plainly intact ORF with no disabling event anywhere near
    # it, and even when both the parent and child of a branch were already
    # individually reported "intact". A tree that is entirely clean below an
    # uncertain root must read as clean, not as one giant uncertain lineage.
    print("root history unresolved by design does not poison confidently intact descendants")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); out = tmp / "o"; out.mkdir()
        pens = tmp / "pens.nwk"; pens.write_text("((A:1,B:1)Node1:1,(C:1,D:1)Node2:1)Root;\n")
        paml = tmp / "paml.nwk"; paml.write_text("((A:1,B:1)P1:1,(C:1,D:1)P2:1)PRoot;\n")
        # Every reconstructed/observed sequence is a clean, identical, intact ORF.
        asr = tmp / "asr.fa"; asr.write_text(">P1\nATGAAACCCGGG\n>P2\nATGAAACCCGGG\n")
        tips = tmp / "tips.fa"; tips.write_text("".join(f">{x}\nATGAAACCCGGG\n" for x in "ABCD"))
        chars = tmp / "chars.tsv"
        write_table(chars, ["character_id", "character_class", "alignment_start", "alignment_end",
                             "length_mod_3", "stop_codon"], [
            {"character_id": "IND0001", "character_class": "indel", "alignment_start": 4,
             "alignment_end": 4, "length_mod_3": 1, "stop_codon": "NA"},
        ])
        # Only the root's state for this frame-shifting-length character is
        # ambiguous (root_history_uncertain=True); every node actually inside
        # the sampled tree has a definite reconstructed state (as a real
        # Sankoff pass always produces), and it is never actually
        # observed/gained on any branch (events.tsv below is empty).
        states = tmp / "states.tsv"
        write_table(states, ["character_id", "node_label", "state"], [
            {"character_id": "IND0001", "node_label": "Root", "state": "ambiguous"},
            {"character_id": "IND0001", "node_label": "Node1", "state": "residue"},
            {"character_id": "IND0001", "node_label": "Node2", "state": "residue"},
        ])
        events = tmp / "events.tsv"
        write_table(events, ["event_id", "branch", "character_class", "biological_interpretation",
                              "length_mod_3", "direction_confident"], [])
        frame = tmp / "frame.tsv"; write_table(frame, ["node_label", "frame_currently_shifted",
                                                         "structural_state_ambiguous"], [])
        r = subprocess.run([sys.executable, str(ROOT / "scripts" / "04_ancestral_orf_walk.py"), "--gene", "T",
                             "--paml-tree", str(paml), "--pensieve-tree", str(pens), "--ancestral-fasta", str(asr),
                             "--tip-alignment", str(tips), "--character-states", str(states), "--characters", str(chars),
                             "--events", str(events), "--frame-arithmetic", str(frame), "--outdir", str(out)],
                            capture_output=True, text=True)
        check(r.returncode == 0, f"uncertain-root synthetic case runs ({r.stderr[:160]})")
        rows = tsv(out / "04_T.orf_transitions_by_branch.tsv")
        br = {x["branch"]: x for x in rows}
        not_intact = [b for b, x in br.items() if x["orf_transition"] != "intact"]
        check(not not_intact,
              f"every branch below the uncertain root resolves to 'intact', not stuck 'uncertain' ({not_intact})")
        check(all(x["child_coding_status"] == "intact" for x in rows),
              "every node's own coding status is intact")
        check(br["Root->Node1"]["transition_evidence"] == "no_inferred_disabling_history"
              and br["Node1->A"]["transition_evidence"] == "no_inferred_disabling_history",
              "resolved branches carry the normal clean-history evidence label, not a lingering uncertainty one")


def test_ambiguous_in_frame_indel_does_not_force_node_uncertain():
    # Real bug found by reviewing genuine HPC runs (GUCA1B/GUCA1C/CNGA3): any
    # character with an ambiguous/unresolved reconstructed state at a node
    # marked that node's whole ORF-completeness call "uncertain", even for a
    # small in-frame indel whose gap/residue polarity cannot possibly change
    # whether the resulting sequence has a valid, complete ORF (it changes
    # neither the reading frame nor the stop-codon content). Genes with many
    # ordinary in-frame indels -- not itself unusual or evidence of anything
    # -- had nearly every internal node marked "uncertain" as a result, which
    # then cascaded into "history_uncertain_no_confirmed_event" everywhere.
    print("an ambiguous in-frame indel does not force a node's ORF status to uncertain")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); out = tmp / "o"; out.mkdir()
        pens = tmp / "pens.nwk"; pens.write_text("((A:1,B:1)Node1:1,C:1)Root;\n")
        paml = tmp / "paml.nwk"; paml.write_text("((A:1,B:1)P1:1,C:1)PRoot;\n")
        asr = tmp / "asr.fa"; asr.write_text(">P1\nATGAAACCCGGG\n")
        tips = tmp / "tips.fa"; tips.write_text("".join(f">{x}\nATGAAACCCGGG\n" for x in "ABC"))
        chars = tmp / "chars.tsv"
        write_table(chars, ["character_id", "character_class", "alignment_start", "alignment_end",
                             "length_mod_3", "stop_codon"], [
            # A 3-bp (in-frame) indel: length_mod_3 = 0.
            {"character_id": "IND0001", "character_class": "indel", "alignment_start": 4,
             "alignment_end": 6, "length_mod_3": 0, "stop_codon": "NA"},
        ])
        # Node1's own reconstruction of this in-frame indel is a genuine
        # parsimony tie (ambiguous) -- e.g. equally good as an insertion on
        # one side or a deletion on the other -- with no bearing on frame or
        # stop-codon content either way.
        states = tmp / "states.tsv"
        write_table(states, ["character_id", "node_label", "state"], [
            {"character_id": "IND0001", "node_label": "Node1", "state": "ambiguous"},
        ])
        events = tmp / "events.tsv"
        write_table(events, ["event_id", "branch", "character_class", "biological_interpretation",
                              "length_mod_3", "direction_confident"], [])
        frame = tmp / "frame.tsv"; write_table(frame, ["node_label", "frame_currently_shifted",
                                                         "structural_state_ambiguous"], [])
        r = subprocess.run([sys.executable, str(ROOT / "scripts" / "04_ancestral_orf_walk.py"), "--gene", "T",
                             "--paml-tree", str(paml), "--pensieve-tree", str(pens), "--ancestral-fasta", str(asr),
                             "--tip-alignment", str(tips), "--character-states", str(states), "--characters", str(chars),
                             "--events", str(events), "--frame-arithmetic", str(frame), "--outdir", str(out)],
                            capture_output=True, text=True)
        check(r.returncode == 0, f"ambiguous in-frame indel synthetic case runs ({r.stderr[:160]})")
        walk = {x["node_label"]: x for x in tsv(out / "04_T.ancestral_orf_walk.tsv")}
        check(walk["Node1"]["coding_status"] == "intact",
              f"Node1 stays intact despite the ambiguous in-frame indel reconstruction ({walk['Node1']})")


def test_sticky_pseudogenic_history():
    print("sticky pseudogenic history after compensatory change")
    # Directly exercise the ORF-walk executable with a deep tree and hand-built
    # event/state tables. Node2->Node1 gets a 1-bp deletion; Node1->A gets a
    # 1-bp residue restoration. A is ORF-intact, but the branch must not turn
    # into a new pseudogenization or erase the inherited history.
    with tempfile.TemporaryDirectory() as td:
        tmp=Path(td); out=tmp/"o"; out.mkdir()
        pens=tmp/"pens.nwk"; pens.write_text("(((A:1,B:1)Node1:1,C:1)Node2:1,D:1)Node3;\n")
        paml=tmp/"paml.nwk"; paml.write_text("(((A:1,B:1)P1:1,C:1)P2:1,D:1)PRoot;\n")
        asr=tmp/"asr.fa"; asr.write_text(">P1\nATGAAACCCGGG\n>P2\nATGAAACCCGGG\n")
        tips=tmp/"tips.fa"; tips.write_text(">A\nATGAAACCCGGG\n>B\nATGAA-CCCGGG\n>C\nATGAAACCCGGG\n>D\nATGAAACCCGGG\n")
        chars=tmp/"chars.tsv"
        write_table(chars,["character_id","character_class","alignment_start","alignment_end","length_mod_3","stop_codon"],[
            {"character_id":"IND0001","character_class":"indel","alignment_start":6,"alignment_end":6,"length_mod_3":1,"stop_codon":"NA"}
        ])
        states=tmp/"states.tsv"
        write_table(states,["character_id","node_label","state"],[
            {"character_id":"IND0001","node_label":"Node3","state":"residue"},
            {"character_id":"IND0001","node_label":"Node2","state":"residue"},
            {"character_id":"IND0001","node_label":"Node1","state":"gap"},
            {"character_id":"IND0001","node_label":"A","state":"residue"},
            {"character_id":"IND0001","node_label":"B","state":"gap"},
            {"character_id":"IND0001","node_label":"C","state":"residue"},
            {"character_id":"IND0001","node_label":"D","state":"residue"},
        ])
        events=tmp/"events.tsv"
        write_table(events,["event_id","branch","character_class","biological_interpretation","length_mod_3","direction_confident"],[
            {"event_id":"del","branch":"Node2->Node1","character_class":"indel","biological_interpretation":"deletion","length_mod_3":1,"direction_confident":True},
            {"event_id":"comp","branch":"Node1->A","character_class":"indel","biological_interpretation":"insertion_or_restoration","length_mod_3":1,"direction_confident":True},
        ])
        frame=tmp/"frame.tsv"; write_table(frame,["node_label","frame_currently_shifted","structural_state_ambiguous"],[])
        r=subprocess.run([sys.executable,str(ROOT/"scripts"/"04_ancestral_orf_walk.py"),"--gene","T",
                          "--paml-tree",str(paml),"--pensieve-tree",str(pens),"--ancestral-fasta",str(asr),
                          "--tip-alignment",str(tips),"--character-states",str(states),"--characters",str(chars),
                          "--events",str(events),"--frame-arithmetic",str(frame),"--outdir",str(out)],
                         capture_output=True,text=True)
        check(r.returncode==0,f"sticky-history synthetic case runs ({r.stderr[:120]})")
        br={x["branch"]:x for x in tsv(out/"04_T.orf_transitions_by_branch.tsv")}
        check(br["Node2->Node1"]["orf_transition"]=="pseudogenization",
              "first confident frameshift is the pseudogenization branch")
        check(br["Node1->A"]["orf_transition"]=="apparent_orf_restoration",
              "compensatory indel is annotated as apparent restoration, not a second loss")
        check(br["Node1->A"]["known_pseudogenic_history"]=="True",
              "compensatory branch retains inherited pseudogenic history")


def test_confident_disabling_events_never_override_a_genuinely_intact_child_orf():
    # Real bug, found by inspecting real PDE6C data directly: Node59 was
    # reported as "pseudogenization" purely because two indel events on its
    # branch were each individually flagged confident/frameshifting (length
    # not a multiple of 3) -- despite Node59's OWN reconstructed sequence,
    # gap-stripped and checked directly, being a genuinely complete, in-frame
    # ORF, consistent with every descendant also being intact. The two
    # events' combined effect on the real sequence was a net multiple of 3
    # (they are not a single clean 3bp event, but together remove exactly
    # 3bp) -- likely reflecting DNA-level rather than strictly codon-level
    # alignment, not a real frameshift ever occurring. The event catalogue is
    # evidence; the actual, directly-checked resulting sequence is ground
    # truth, and must win.
    print("catalogued confident-disabling events never override a child whose own reconstructed ORF is genuinely complete")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); out = tmp / "o"; out.mkdir()
        pens = tmp / "pens.nwk"
        pens.write_text("(((ChildTip1:1,ChildTip2:1)Child:1,Sibling:1)Parent:1,Outgroup:1)Root;\n")
        paml = tmp / "paml.nwk"
        paml.write_text("(((ChildTip1:1,ChildTip2:1)P2:1,Sibling:1)P1:1,Outgroup:1)PRoot;\n")
        # Parent's own real, complete 15bp ORF: ATG-CAA-CCA-CCC-TGA.
        # Child's PAML scaffold is the same raw sequence; two indel
        # characters (1bp at column 4, 2bp at columns 7-8 -- each on its own
        # not a multiple of 3) get gapped out for Child specifically,
        # leaving exactly ATG-AAA-CCC-TGA: a genuinely complete, in-frame
        # 12bp ORF once the 3 gapped columns are stripped.
        asr = tmp / "asr.fa"
        asr.write_text(">P1\nATGCAACCACCCTGA\n>P2\nATGCAACCACCCTGA\n")
        # Tip alignment width (15) must match the PAML scaffold width, exactly
        # as it always does in a real canonical-alignment run -- ChildTip1/2's
        # own real sequence is gapped at the same 3 columns removed for Child
        # above (so their gap-stripped real ending is still ATGAAACCCTGA);
        # Sibling/Outgroup keep the full, un-gapped 15bp sequence.
        tips = tmp / "tips.fa"
        tips.write_text(
            ">ChildTip1\nATG-AA--ACCCTGA\n>ChildTip2\nATG-AA--ACCCTGA\n"
            ">Sibling\nATGCAACCACCCTGA\n>Outgroup\nATGCAACCACCCTGA\n"
        )
        chars = tmp / "chars.tsv"
        write_table(chars, ["character_id", "character_class", "alignment_start", "alignment_end",
                             "length_mod_3", "stop_codon"], [
            {"character_id": "IND0001", "character_class": "indel", "alignment_start": 4,
             "alignment_end": 4, "length_mod_3": 1, "stop_codon": "NA"},
            {"character_id": "IND0002", "character_class": "indel", "alignment_start": 7,
             "alignment_end": 8, "length_mod_3": 2, "stop_codon": "NA"},
        ])
        states = tmp / "states.tsv"
        write_table(states, ["character_id", "node_label", "state"], [
            {"character_id": "IND0001", "node_label": "Root", "state": "residue"},
            {"character_id": "IND0001", "node_label": "Parent", "state": "residue"},
            {"character_id": "IND0001", "node_label": "Child", "state": "gap"},
            {"character_id": "IND0002", "node_label": "Root", "state": "residue"},
            {"character_id": "IND0002", "node_label": "Parent", "state": "residue"},
            {"character_id": "IND0002", "node_label": "Child", "state": "gap"},
        ])
        events = tmp / "events.tsv"
        write_table(events, ["event_id", "branch", "character_class", "biological_interpretation",
                              "length_mod_3", "direction_confident"], [
            {"event_id": "ev1", "branch": "Parent->Child", "character_class": "indel",
             "biological_interpretation": "deletion", "length_mod_3": 1, "direction_confident": True},
            {"event_id": "ev2", "branch": "Parent->Child", "character_class": "indel",
             "biological_interpretation": "deletion", "length_mod_3": 2, "direction_confident": True},
        ])
        frame = tmp / "frame.tsv"
        write_table(frame, ["node_label", "frame_currently_shifted", "structural_state_ambiguous"], [])
        r = subprocess.run([sys.executable, str(ROOT / "scripts" / "04_ancestral_orf_walk.py"), "--gene", "T",
                             "--paml-tree", str(paml), "--pensieve-tree", str(pens), "--ancestral-fasta", str(asr),
                             "--tip-alignment", str(tips), "--character-states", str(states), "--characters", str(chars),
                             "--events", str(events), "--frame-arithmetic", str(frame), "--outdir", str(out)],
                            capture_output=True, text=True)
        check(r.returncode == 0, f"compensating-frameshift synthetic case runs ({r.stderr[:160]})")
        walk = {r["node_label"]: r for r in tsv(out / "04_T.ancestral_orf_walk.tsv")}
        check(walk.get("Child", {}).get("coding_status") == "intact",
              f"Child's own reconstructed sequence is correctly a complete ORF ({walk.get('Child')})")
        br = {x["branch"]: x for x in tsv(out / "04_T.orf_transitions_by_branch.tsv")}
        check(br["Parent->Child"]["orf_transition"] == "intact",
              f"the branch is intact despite two catalogued confident-disabling events, "
              f"because Child's own ORF is genuinely complete ({br['Parent->Child']})")
        check(br["Parent->Child"]["confident_disabling_events_on_branch"] == "ev1,ev2",
              "the catalogued events remain visible in the audit trail even though they didn't win")


def test_sparse_root_stop_evidence_does_not_poison_whole_tree():
    # Real bug, found by inspecting real GUCY2F and CNGB3 data directly: a
    # stop_mask character observed "present" in as few as 1-3 tips (out of
    # 100+), with n_observed_absent == 0 everywhere else, still gets a Fitch
    # root reconstruction of "stop_present" -- with zero contradicting
    # "absent" evidence anywhere, "present" costs zero changes at every node
    # including the root, so parsimony reports a perfect score even though
    # only a handful of, likely deeply nested, tips actually carry it. Before
    # the fix, that single unsupported character flipped root_disabling True
    # for the WHOLE gene, so every independent pseudogenization anywhere in
    # the tree rendered as "already_pseudogenic" (inherited) instead of its
    # own fresh "pseudogenization" call -- both real genes had zero
    # "pseudogenization" branches and only pale/inherited ones. A stop
    # character only has real phylogenetic contrast, and can only genuinely
    # support a pre-root loss, if at least one tip is confidently observed
    # WITHOUT it (n_observed_absent > 0).
    print("sparse root stop_mask evidence (n_observed_absent == 0) must not force root_disabling")

    def run_case(n_observed_absent, extra_state):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td); out = tmp / "o"; out.mkdir()
            pens = tmp / "pens.nwk"
            pens.write_text("(((A:1,B:1)Node1:1,C:1)Node2:1,D:1)Node3;\n")
            paml = tmp / "paml.nwk"
            paml.write_text("(((A:1,B:1)P1:1,C:1)P2:1,D:1)PRoot;\n")
            asr = tmp / "asr.fa"
            asr.write_text(">P1\nATGAAACCCGGG\n>P2\nATGAAACCCGGG\n")
            tips = tmp / "tips.fa"
            tips.write_text(">A\nATGAA-CCCGGG\n>B\nATGAAACCCGGG\n>C\nATGAAACCCGGG\n>D\nATGAAACCCGGG\n")
            chars = tmp / "chars.tsv"
            write_table(chars, ["character_id", "character_class", "alignment_start", "alignment_end",
                                 "length_mod_3", "stop_codon", "n_observed_present", "n_observed_absent"], [
                {"character_id": "IND0001", "character_class": "indel", "alignment_start": 6,
                 "alignment_end": 6, "length_mod_3": 1, "stop_codon": "NA",
                 "n_observed_present": "NA", "n_observed_absent": "NA"},
                {"character_id": "STOPQ", "character_class": "stop_mask", "alignment_start": 20,
                 "alignment_end": 22, "length_mod_3": 0, "stop_codon": "TAA",
                 "n_observed_present": 1, "n_observed_absent": n_observed_absent},
            ])
            states = tmp / "states.tsv"
            rows = [
                {"character_id": "IND0001", "node_label": "Node3", "state": "residue"},
                {"character_id": "IND0001", "node_label": "Node2", "state": "residue"},
                {"character_id": "IND0001", "node_label": "Node1", "state": "residue"},
                {"character_id": "IND0001", "node_label": "A", "state": "gap"},
                {"character_id": "IND0001", "node_label": "B", "state": "residue"},
                {"character_id": "IND0001", "node_label": "C", "state": "residue"},
                {"character_id": "IND0001", "node_label": "D", "state": "residue"},
                {"character_id": "STOPQ", "node_label": "Node3", "state": "stop_present"},
            ]
            if extra_state:
                rows.append(extra_state)
            write_table(states, ["character_id", "node_label", "state"], rows)
            events = tmp / "events.tsv"
            write_table(events, ["event_id", "branch", "character_class", "biological_interpretation",
                                  "length_mod_3", "direction_confident"], [
                {"event_id": "del", "branch": "Node1->A", "character_class": "indel",
                 "biological_interpretation": "deletion", "length_mod_3": 1, "direction_confident": True},
            ])
            frame = tmp / "frame.tsv"
            write_table(frame, ["node_label", "frame_currently_shifted", "structural_state_ambiguous"], [])
            r = subprocess.run([sys.executable, str(ROOT / "scripts" / "04_ancestral_orf_walk.py"), "--gene", "T",
                                 "--paml-tree", str(paml), "--pensieve-tree", str(pens), "--ancestral-fasta", str(asr),
                                 "--tip-alignment", str(tips), "--character-states", str(states), "--characters", str(chars),
                                 "--events", str(events), "--frame-arithmetic", str(frame), "--outdir", str(out)],
                                capture_output=True, text=True)
            check(r.returncode == 0, f"root-stop-evidence synthetic case runs ({r.stderr[:160]})")
            br = {x["branch"]: x for x in tsv(out / "04_T.orf_transitions_by_branch.tsv")}
            return br["Node1->A"]["orf_transition"]

    # Both "pseudogenization" (parent history confidently clean) and
    # "confirmed_disabling_event_first_loss_unresolved" (parent history not
    # confidently resolved either way, but this branch's own event still
    # stands on its own) render identically as solid red in 05_plot_events.R
    # -- what must NOT happen is the sticky pale "already_pseudogenic".
    sparse = run_case(0, None)
    check(sparse in ("pseudogenization", "confirmed_disabling_event_first_loss_unresolved"),
          f"a stop_mask character present in 1 tip with zero confidently-absent tips does not make an "
          f"unrelated, independently-disabled branch look 'already_pseudogenic' (got {sparse!r})")

    supported = run_case(1, {"character_id": "STOPQ", "node_label": "D", "state": "stop_absent"})
    check(supported == "already_pseudogenic",
          f"a stop_mask character with a genuinely confirmed-absent tip still correctly forces sticky "
          f"root-level pseudogenic history (got {supported!r})")


def test_intact_parent_resets_history_for_independently_disabled_children():
    # Real bug, found by inspecting real GUCY2F/CNGB3 data: sticky
    # known_pseudogenic_history (see test_sticky_pseudogenic_history) is
    # deliberately True for a node whose OWN current sequence is intact but
    # whose more distant ancestor was pseudogenized ("apparent_orf_restoration").
    # That stickiness is correct for THAT node's own record, but blindly
    # propagating it to its children made two sister branches, each with
    # their own brand-new independent disabling event off a genuinely intact
    # common ancestor, render as "already_pseudogenic" (pale, inherited)
    # instead of their own fresh "pseudogenization" (solid red) call --
    # visually contradicting the grey/intact parent right above them.
    print("children of a currently-intact parent get fresh pseudogenization calls, not inherited sticky history")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); out = tmp / "o"; out.mkdir()
        pens = tmp / "pens.nwk"
        pens.write_text("((((A2a:1,A2b:1)A:1,B:1)Node1:1,C:1)Node2:1,D:1)Node3;\n")
        paml = tmp / "paml.nwk"
        paml.write_text("((((A2a:1,A2b:1)P1a:1,B:1)P1:1,C:1)P2:1,D:1)PRoot;\n")
        asr = tmp / "asr.fa"
        asr.write_text(">P1a\nATGAAACCCGGG\n>P1\nATGAAACCCGGG\n>P2\nATGAAACCCGGG\n")
        tips = tmp / "tips.fa"
        tips.write_text(
            ">A2a\nATGAAACC-GGG\n>A2b\nATGAAACCCG-G\n>B\nATGAA-CCCGGG\n"
            ">C\nATGAAACCCGGG\n>D\nATGAAACCCGGG\n"
        )
        chars = tmp / "chars.tsv"
        write_table(chars, ["character_id", "character_class", "alignment_start", "alignment_end",
                             "length_mod_3", "stop_codon"], [
            {"character_id": "IND0001", "character_class": "indel", "alignment_start": 6,
             "alignment_end": 6, "length_mod_3": 1, "stop_codon": "NA"},
            {"character_id": "IND0002", "character_class": "indel", "alignment_start": 9,
             "alignment_end": 9, "length_mod_3": 1, "stop_codon": "NA"},
            {"character_id": "IND0003", "character_class": "indel", "alignment_start": 11,
             "alignment_end": 11, "length_mod_3": 1, "stop_codon": "NA"},
        ])
        states = tmp / "states.tsv"
        write_table(states, ["character_id", "node_label", "state"], [
            # IND0001: Node2->Node1 deletion, Node1->A compensating restoration
            # (the same shape as test_sticky_pseudogenic_history).
            {"character_id": "IND0001", "node_label": "Node3", "state": "residue"},
            {"character_id": "IND0001", "node_label": "Node2", "state": "residue"},
            {"character_id": "IND0001", "node_label": "Node1", "state": "gap"},
            {"character_id": "IND0001", "node_label": "A", "state": "residue"},
            {"character_id": "IND0001", "node_label": "A2a", "state": "residue"},
            {"character_id": "IND0001", "node_label": "A2b", "state": "residue"},
            {"character_id": "IND0001", "node_label": "B", "state": "gap"},
            {"character_id": "IND0001", "node_label": "C", "state": "residue"},
            {"character_id": "IND0001", "node_label": "D", "state": "residue"},
            # IND0002: A2a's own new, independent deletion.
            {"character_id": "IND0002", "node_label": "A", "state": "residue"},
            {"character_id": "IND0002", "node_label": "A2a", "state": "gap"},
            {"character_id": "IND0002", "node_label": "A2b", "state": "residue"},
            # IND0003: A2b's own new, independent deletion.
            {"character_id": "IND0003", "node_label": "A", "state": "residue"},
            {"character_id": "IND0003", "node_label": "A2a", "state": "residue"},
            {"character_id": "IND0003", "node_label": "A2b", "state": "gap"},
        ])
        events = tmp / "events.tsv"
        write_table(events, ["event_id", "branch", "character_class", "biological_interpretation",
                              "length_mod_3", "direction_confident"], [
            {"event_id": "del", "branch": "Node2->Node1", "character_class": "indel",
             "biological_interpretation": "deletion", "length_mod_3": 1, "direction_confident": True},
            {"event_id": "comp", "branch": "Node1->A", "character_class": "indel",
             "biological_interpretation": "insertion_or_restoration", "length_mod_3": 1, "direction_confident": True},
            {"event_id": "newA2a", "branch": "A->A2a", "character_class": "indel",
             "biological_interpretation": "deletion", "length_mod_3": 1, "direction_confident": True},
            {"event_id": "newA2b", "branch": "A->A2b", "character_class": "indel",
             "biological_interpretation": "deletion", "length_mod_3": 1, "direction_confident": True},
        ])
        frame = tmp / "frame.tsv"
        write_table(frame, ["node_label", "frame_currently_shifted", "structural_state_ambiguous"], [])
        r = subprocess.run([sys.executable, str(ROOT / "scripts" / "04_ancestral_orf_walk.py"), "--gene", "T",
                             "--paml-tree", str(paml), "--pensieve-tree", str(pens), "--ancestral-fasta", str(asr),
                             "--tip-alignment", str(tips), "--character-states", str(states), "--characters", str(chars),
                             "--events", str(events), "--frame-arithmetic", str(frame), "--outdir", str(out)],
                            capture_output=True, text=True)
        check(r.returncode == 0, f"intact-parent synthetic case runs ({r.stderr[:200]})")
        walk = {r["node_label"]: r for r in tsv(out / "04_T.ancestral_orf_walk.tsv")}
        check(walk.get("A", {}).get("coding_status") == "intact",
              f"A's own reconstructed sequence is intact despite inherited sticky history ({walk.get('A')})")
        br = {x["branch"]: x for x in tsv(out / "04_T.orf_transitions_by_branch.tsv")}
        check(br["Node1->A"]["orf_transition"] == "apparent_orf_restoration",
              "A itself is still correctly labelled apparent_orf_restoration")
        check(br["A->A2a"]["orf_transition"] == "pseudogenization",
              f"A2a gets its own fresh pseudogenization call, not inherited history ({br.get('A->A2a')})")
        check(br["A->A2b"]["orf_transition"] == "pseudogenization",
              f"A2b gets its own fresh pseudogenization call, not inherited history ({br.get('A->A2b')})")


def test_map_raw_nt_to_aln_skips_macse_synthetic_ambiguity_code():
    # Real bug, found on real data (GUCA1C/Desmodus_rotundus, reported
    # directly by the user). map_raw_nt_to_aln() used to treat every
    # non-gap, non-'!' character in a MACSE alignment row as a real raw
    # nucleotide. MACSE's own output can contain a literal ambiguity code
    # like 'N' as a synthetic insertion -- not present in the raw input at
    # all -- with no '!' involved. Treating it as real silently miscounted
    # by one: in the real case, this shifted the 30-codon conserved-block
    # boundary one raw nucleotide early for Desmodus_rotundus, so its real
    # nucleotide exactly at the boundary (raw position 90, a 'C') fell
    # between the accepted block and the remainder and was silently
    # dropped -- then reported as a false single-nucleotide deletion event
    # unique to that species.
    print("map_raw_nt_to_aln correctly excludes a MACSE-synthetic 'N', not just '!'")
    mod = load("01_run_macse_and_extract_events")
    raw_seq = "ATGAAACCCGGG"  # 12 real raw nucleotides, no ambiguity
    # A synthetic leading 'N' (matching no raw character) plus a real
    # internal gap -- exactly the real-data shape (MACSE row for Desmodus
    # started 'NTGGGC------AAT---...').
    aligned = "N" + "ATGAAA" + "-" + "CCCGGG"
    col = mod.map_raw_nt_to_aln(aligned, 6, raw_seq)
    check(col == 7, f"synthetic leading 'N' is excluded from the raw-nucleotide count (got column {col}, expected 7)")
    col_last = mod.map_raw_nt_to_aln(aligned, 12, raw_seq)
    check(col_last == len(aligned),
          f"the last real raw nucleotide still maps to the true final alignment column (got {col_last})")

    # Control: a genuine 'N' actually present IN the raw sequence itself
    # (real ambiguity in the original data, not a MACSE insertion) must
    # still be counted as a real raw nucleotide, not skipped.
    raw_with_n = "ATGNAACCCGGG"
    aligned_with_n = "N" + "ATGNAA" + "-" + "CCCGGG"
    col_genuine = mod.map_raw_nt_to_aln(aligned_with_n, 4, raw_with_n)
    check(col_genuine == 5,
          f"a genuine 'N' present in the raw sequence itself is still counted, not treated as synthetic (got {col_genuine})")


def test_raw_to_alignment_map_survives_multiple_scattered_insertions():
    # Real bug, found running the v4.3 fix above on real GUCA1C data: a
    # naive single-pass greedy character match (advance on match, else
    # skip) correctly handles ONE synthetic MACSE insertion, but perma-
    # nently desynchronises as soon as there is a SECOND, independent one
    # later in the same row -- once a real character fails to match by
    # coincidence, greedy matching has no way to recover, and can stall
    # with zero further matches for the rest of the sequence. On real data
    # this blocked GUCA1C entirely: Desmodus_rotundus's own row had a
    # synthetic insertion near the very start (the conserved-block
    # boundary case) AND a second, unrelated one further into the gene,
    # and greedy matching found no STOP-coordinate mapping anywhere past
    # the first one. build_raw_to_alignment_map() uses a proper sequence
    # alignment (difflib) instead, which does not have this failure mode.
    print("raw-to-alignment map survives more than one synthetic insertion in the same row")
    mod = load("01_run_macse_and_extract_events")
    # A non-repetitive 30nt sequence (fixed seed) -- a sequence built from
    # short repeated motifs (e.g. homopolymer-like runs) can give a generic
    # sequence matcher more than one equally-short valid alignment to choose
    # between, which is a property of the test data, not of real biological
    # sequences at this length; a realistic-looking sequence avoids that
    # confound entirely.
    raw_seq = "AAGTTTAACAAATTATCTCGATCGGTTGAA"  # len 30, fixed-seed random
    # Two independent synthetic insertions: an 'N' 10nt in, and a '!'
    # partial-codon marker another 10nt later -- exactly the real
    # GUCA1C/Desmodus_rotundus shape (an insertion at the conserved-block
    # boundary, plus one more, unrelated, later in the gene).
    aligned = raw_seq[:10] + "N" + raw_seq[10:20] + "!" + raw_seq[20:]
    raw_map = mod.build_raw_to_alignment_map(aligned, raw_seq)
    check(len(raw_map) == len(raw_seq),
          f"every real raw position is mapped despite two separate synthetic insertions "
          f"(mapped {len(raw_map)} of {len(raw_seq)})")
    # The last real raw nucleotide must still be found, not lost to a stall
    # after the first insertion.
    check(len(raw_seq) in raw_map,
          f"the last real raw nucleotide is still mapped (raw_map keys: {sorted(raw_map)})")
    # Spot-check that mapped columns actually point at the matching character.
    for raw_pos, col in raw_map.items():
        check(aligned[col - 1].upper() == raw_seq[raw_pos - 1].upper(),
              f"raw position {raw_pos} maps to column {col}, whose character "
              f"({aligned[col-1]!r}) matches the raw base ({raw_seq[raw_pos-1]!r})")


def test_no_muscle_in_active_runtime_code():
    # Mandatory acceptance criterion (v4.2): MUSCLE is removed completely.
    # There must be no --aligner option and no runtime code path that can
    # invoke an executable literally named "muscle". Historical CHANGELOG
    # prose describing the removed feature is fine; a live subprocess call
    # or CLI flag is not.
    print("no active runtime code invokes muscle or exposes --aligner")
    active_dirs = [ROOT / "bin", ROOT / "scripts", ROOT / "templates"]
    hits = []
    for d in active_dirs:
        if not d.exists():
            continue
        for path in d.rglob("*"):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            try:
                text = path.read_text()
            except (UnicodeDecodeError, OSError):
                continue
            for lineno, line in enumerate(text.splitlines(), start=1):
                if "muscle" in line.lower() or "--aligner" in line:
                    hits.append(f"{path.relative_to(ROOT)}:{lineno}: {line.strip()}")
    check(not hits, f"no muscle/--aligner reference in bin/scripts/templates; found: {hits[:10]}")


def main():
    test_static_and_cli()
    test_terminal_stop_stripped_from_every_input_sequence()
    test_breakpoint_decomposition()
    test_long_deletion_does_not_manufacture_a_spurious_shared_core_event()
    test_functional_consensus_fixes_shared_indel_as_ancestral()
    test_guca1b_breakpoint_history_end_to_end()
    test_contiguous_same_type_indel_fragments_are_merged_into_one_event()
    test_adjacent_different_type_indel_events_are_not_merged()
    test_contiguous_shared_events_with_identical_affected_tips_are_merged()
    test_events_with_different_affected_tips_are_not_merged_even_on_the_same_branch()
    test_overlapping_affected_tips_contiguous_events_are_clubbed()
    test_nested_interior_event_does_not_block_merging_the_events_around_it()
    test_parsimony_tie_and_direction()
    test_stop_alleles_are_separate()
    test_shared_stop_not_missed_when_only_one_tip_is_registered()
    test_compensated_frameshift_stop_classification()
    test_build_stop_registry_rejects_non_contiguous_codon_mapping()
    test_build_stop_registry_rejects_off_frame_codon_coincidence()
    test_frame_dependent_stop_not_resurrected_by_coincidental_allele_match()
    test_canonical_alignment_prepare()
    test_defined_alignment_direct_codon_stop_scan_without_macse()
    test_prepare_defined_alignment_frame_fix_and_masking()
    test_phylogenetic_sequence_ordering()
    test_runner_orchestration_through_events()
    test_runner_orchestration_through_integrate()
    test_mock_paml_integration_and_root_policy()
    test_partial_orf_status_distinguishes_truncation_from_pseudogenization()
    test_root_adjacent_uncertainty_is_explicit()
    test_uncertain_root_resolves_to_intact_once_descendants_are_confidently_intact()
    test_ambiguous_in_frame_indel_does_not_force_node_uncertain()
    test_sticky_pseudogenic_history()
    test_confident_disabling_events_never_override_a_genuinely_intact_child_orf()
    test_sparse_root_stop_evidence_does_not_poison_whole_tree()
    test_intact_parent_resets_history_for_independently_disabled_children()
    test_map_raw_nt_to_aln_skips_macse_synthetic_ambiguity_code()
    test_raw_to_alignment_map_survives_multiple_scattered_insertions()
    test_no_muscle_in_active_runtime_code()
    if FAILURES:
        print(f"\n{len(FAILURES)} failure(s)")
        sys.exit(1)
    print("\nBackend consistency test passed.")


if __name__ == "__main__":
    main()
