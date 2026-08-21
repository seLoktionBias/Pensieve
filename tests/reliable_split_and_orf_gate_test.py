#!/usr/bin/env python3
"""ORF-based MACSE input split + the pre-PAML complete-ORF gate (v4.8).

Replaces the conserved-start-block heuristic, which is removed. That heuristic
detected MACSE corrupting an intact reading frame only AFTER the fact, only in
the first 30 codons, only for genes with >=5 complete ORFs, and only at the 5'
end. Instead, step 00 now splits the gapless CDS by ORF completeness and step 01
tells MACSE which sequences it may trust:

    macse -prog alignSequences -seq complete_seqs.fa -seq_lr incomplete_seqs.fa \\
          -fs 1000 -fs_term 1000

so a frameshift proposed inside a complete ORF is prohibitively expensive
everywhere in the gene, while frameshifts in the less-reliable set -- where the
real lesions are -- keep MACSE's default cost.

Because that is a claim about the output, step 02 enforces it as a hard gate
before anything reaches PAML: every complete ORF must reproduce its input
exactly after degapping, carry zero MACSE '!' markers, and occupy only whole
codon cells (three residues or '---').
"""
from __future__ import annotations

import csv
import importlib.util
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


def load(script):
    spec = importlib.util.spec_from_file_location(Path(script).stem.strip("0123456789_"),
                                                  ROOT / "scripts" / script)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fasta(seqs):
    return "".join(f">{k}\n{v}\n" for k, v in seqs.items())


def make_mock_macse(path, argv_log):
    """A stand-in for MACSE that records its argv and emits a padded alignment."""
    path.write_text(f"""#!/usr/bin/env python3
import sys
open({str(argv_log)!r}, "w").write("\\n".join(sys.argv[1:]))
args = sys.argv[1:]
def val(k): return args[args.index(k) + 1]
seqs = {{}}
for flag in ("-seq", "-seq_lr"):
    if flag in args:
        name = None
        for line in open(val(flag)):
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]; seqs[name] = ""
            elif name:
                seqs[name] += line
width = max(len(v) for v in seqs.values())
width += (-width) % 3
with open(val("-out_NT"), "w") as out:
    for k, v in seqs.items():
        out.write(">%s\\n%s\\n" % (k, v + "-" * (width - len(v))))
with open(val("-out_AA"), "w") as out:
    for k in seqs:
        out.write(">%s\\n%s\\n" % (k, "X" * (width // 3)))
""")
    path.chmod(0o755)
    return path


# ---------------------------------------------------------------- removal
print("the conserved-start-block heuristic is gone from the runtime code")
banned = ("detect_conserved_block", "resolve_conserved_block", "realign_block_content",
          "block_reference", "block_quality_scores", "conserved_block_correction")
hits = []
for d in (ROOT / "bin", ROOT / "scripts", ROOT / "templates"):
    for path in d.rglob("*"):
        if not path.is_file() or "__pycache__" in path.parts:
            continue
        try:
            text = path.read_text()
        except (UnicodeDecodeError, OSError):
            continue
        for name in banned:
            if name in text:
                hits.append(f"{path.relative_to(ROOT)}:{name}")
check(not hits, f"no conserved-block symbol survives in bin/scripts/templates; found: {hits}")

# ---------------------------------------------------------------- step 00 split
print("step 00 splits the gapless CDS by ORF completeness")
COMPLETE = "ATGAAACCCGGGTTTAAACCCGGGTAA"          # ATG, len%3==0, no internal stop
BROKEN   = "ATGAAACCCTAAGGGTTTAAACCCGGGTAA"       # internal in-frame stop
NO_ATG   = "GGGAAACCCGGGTTTAAACCCGGGTAA"          # no start codon
with tempfile.TemporaryDirectory() as d:
    tmp = Path(d)
    seqs = {"A": COMPLETE, "B": COMPLETE, "C": BROKEN, "D": NO_ATG}
    (tmp / "in.fa").write_text(fasta(seqs))
    (tmp / "in.nwk").write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "00_prune_and_check_orf.py"),
                        "--gene", "T", "--fasta", str(tmp / "in.fa"), "--tree", str(tmp / "in.nwk"),
                        "--outdir", str(tmp / "out")], capture_output=True, text=True)
    check(r.returncode == 0, f"step 00 runs ({r.stderr[-300:]})")
    comp = tmp / "out" / "00_T.complete_seqs.fa"
    inc = tmp / "out" / "00_T.incomplete_seqs.fa"
    check(comp.exists() and inc.exists(), "both complete_seqs.fa and incomplete_seqs.fa are written")
    if comp.exists() and inc.exists():
        cnames = [l[1:].strip() for l in comp.read_text().splitlines() if l.startswith(">")]
        inames = [l[1:].strip() for l in inc.read_text().splitlines() if l.startswith(">")]
        check(sorted(cnames) == ["A", "B"], f"complete ORFs go to complete_seqs.fa (got {sorted(cnames)})")
        check(sorted(inames) == ["C", "D"],
              f"an internal stop and a missing ATG both go to incomplete_seqs.fa (got {sorted(inames)})")
        status = {r_["species"]: r_ for r_ in tsv(tmp / "out" / "00_T.orf_status.tsv")}
        check(all((status[n]["complete_orf"] == "True") == (n in cnames) for n in status),
              "the split agrees with orf_status.tsv's own complete_orf column")
        body = "".join(l for l in comp.read_text().splitlines() if not l.startswith(">"))
        check("-" not in body, "the split files are gapless (MACSE receives unaligned CDS)")
        # The terminal stop must already be gone -- step 00 strips it before the split.
        check(not body.startswith(COMPLETE),
              "the terminal stop is stripped before the split, not carried into MACSE input")

# ---------------------------------------------------------------- step 01 command
print("step 01 builds the reliable/less-reliable MACSE command")
with tempfile.TemporaryDirectory() as d:
    tmp = Path(d)
    s00 = tmp / "r00"; s00.mkdir()
    s01 = tmp / "r01"; s01.mkdir()
    (s00 / "00_T.complete_seqs.fa").write_text(fasta({"A": COMPLETE, "B": COMPLETE}))
    (s00 / "00_T.incomplete_seqs.fa").write_text(fasta({"C": BROKEN}))
    (s00 / "00_T.common_species.gapless_for_macse.fasta").write_text(
        fasta({"A": COMPLETE, "B": COMPLETE, "C": BROKEN}))
    (s00 / "00_T.common_species.tree").write_text("((A:1,B:1),C:1);\n")
    with open(s00 / "00_T.orf_status.tsv", "w") as h:
        h.write("gene\tspecies\tcomplete_orf\n")
        for n, v in (("A", True), ("B", True), ("C", False)):
            h.write(f"T\t{n}\t{v}\n")
    (s00 / "00_T.orf_failures.tsv").write_text("gene\tspecies\tfailure_type\n")
    argv_log = tmp / "argv.txt"
    mock = make_mock_macse(tmp / "mock_macse", argv_log)
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "01_run_macse_and_extract_events.py"),
                        "--gene", "T", "--step00-dir", str(s00), "--outdir", str(s01),
                        "--macse-cmd", str(mock)], capture_output=True, text=True)
    check(r.returncode == 0, f"step 01 runs ({r.stderr[-400:]})")
    argv = argv_log.read_text().splitlines() if argv_log.exists() else []
    check("-seq" in argv and argv[argv.index("-seq") + 1].endswith("00_T.complete_seqs.fa"),
          f"-seq is the complete-ORF file (argv={argv})")
    check("-seq_lr" in argv and argv[argv.index("-seq_lr") + 1].endswith("00_T.incomplete_seqs.fa"),
          "-seq_lr is the incomplete-ORF file")
    check("-fs" in argv and argv[argv.index("-fs") + 1] == "1000",
          f"-fs 1000 is passed (got {argv[argv.index('-fs') + 1] if '-fs' in argv else 'absent'})")
    check("-fs_term" in argv and argv[argv.index("-fs_term") + 1] == "1000",
          "-fs_term 1000 is passed")
    check(not (s01 / "01_T.conserved_block_correction.tsv").exists(),
          "no conserved-block correction table is written any more")

print("a gene with no complete ORF falls back to the undivided input")
with tempfile.TemporaryDirectory() as d:
    tmp = Path(d)
    s00 = tmp / "r00"; s00.mkdir()
    s01 = tmp / "r01"; s01.mkdir()
    (s00 / "00_T.complete_seqs.fa").write_text("")
    (s00 / "00_T.incomplete_seqs.fa").write_text(fasta({"C": BROKEN, "D": NO_ATG}))
    (s00 / "00_T.common_species.gapless_for_macse.fasta").write_text(fasta({"C": BROKEN, "D": NO_ATG}))
    (s00 / "00_T.common_species.tree").write_text("(C:1,D:1);\n")
    with open(s00 / "00_T.orf_status.tsv", "w") as h:
        h.write("gene\tspecies\tcomplete_orf\n"); h.write("T\tC\tFalse\nT\tD\tFalse\n")
    (s00 / "00_T.orf_failures.tsv").write_text("gene\tspecies\tfailure_type\n")
    argv_log = tmp / "argv.txt"
    mock2 = make_mock_macse(tmp / "mock_macse", argv_log)
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "01_run_macse_and_extract_events.py"),
                        "--gene", "T", "--step00-dir", str(s00), "--outdir", str(s01),
                        "--macse-cmd", str(mock2)], capture_output=True, text=True)
    check(r.returncode == 0, f"step 01 runs with no complete ORF ({r.stderr[-300:]})")
    argv = argv_log.read_text().splitlines() if argv_log.exists() else []
    check("-seq_lr" not in argv, f"-seq_lr alone is never passed to MACSE (argv={argv})")
    check("-fs" not in argv,
          "frameshifts are NOT forbidden when every sequence is known to be frameshifted")

# ---------------------------------------------------------------- step 02 gate
print("the pre-PAML gate accepts an intact complete ORF and rejects a corrupted one")
prep = load("02_prepare_asr_inputs.py")
with tempfile.TemporaryDirectory() as d:
    out = Path(d)
    inputs = {"A": "ATGAAACCCGGG", "B": "ATGAAACCC", "C": "ATGAAACCCGGG"}
    order = ["A", "B", "C"]

    good_native = {"A": "ATGAAACCCGGG", "B": "ATGAAACCC---", "C": "ATGAAACCCGGG"}
    rows = prep.validate_complete_orf_alignment(
        "T", dict(good_native), good_native, order, ["A", "B"], inputs, out, "test")
    verdicts = {r["species"]: r["verdict"] for r in rows}
    check(verdicts["A"] == "pass" and verdicts["B"] == "pass",
          f"an intact complete ORF passes (got {verdicts})")
    check(verdicts["C"] == "not_applicable",
          "an incomplete-ORF sequence is reported but never gates the pipeline")

    def expect_fail(label, native, source, reason_fragment):
        try:
            prep.validate_complete_orf_alignment(
                "T", source, native, order, ["A", "B"], inputs, out, "test")
        except SystemExit as exc:
            report = out / "02_T.complete_orf_alignment_validation.tsv"
            reasons = " ".join(r["failure_reason"] for r in tsv(report))
            check(reason_fragment in reasons,
                  f"{label}: the recorded reason names the violation (got {reasons!r})")
            check(str(report) in str(exc), f"{label}: the error points at the diagnostic report")
            return
        check(False, f"{label}: should have stopped the pipeline, but did not")

    # 1. a real nucleotide altered by the aligner
    bad = dict(good_native); bad["A"] = "ATGAAANCCGGG"
    expect_fail("altered residue", bad, dict(bad), "degapped_row_does_not_reproduce_input")

    # 2. a MACSE frameshift marker inside a complete ORF
    src = dict(good_native); src["A"] = "ATGAAACCCGG!"
    nat = {k: v.replace("!", "-") for k, v in src.items()}
    inputs2 = dict(inputs); inputs2["A"] = "ATGAAACCCGG"
    try:
        prep.validate_complete_orf_alignment("T", src, nat, order, ["A", "B"], inputs2, out, "test")
        check(False, "frameshift marker: should have stopped the pipeline, but did not")
    except SystemExit:
        reasons = " ".join(r["failure_reason"] for r in tsv(out / "02_T.complete_orf_alignment_validation.tsv"))
        check("frameshift_marker" in reasons, f"frameshift marker is named as the violation ({reasons!r})")

    # 3. a codon split across a gap boundary
    bad = {"A": "ATGAAACC-CGGG", "B": "ATGAAACCC----", "C": "ATGAAACCCGGG-"}
    inputs3 = {"A": "ATGAAACCCGGG", "B": "ATGAAACCC", "C": "ATGAAACCCGGG"}
    expect_fail("split codon", bad, dict(bad), "incomplete_codon_cell")

    # The report survives a failure -- that is the point of it.
    check((out / "02_T.complete_orf_alignment_validation.tsv").exists(),
          "the diagnostic report is on disk after a failure")

if FAILURES:
    print(f"\n{len(FAILURES)} failure(s):")
    for f in FAILURES:
        print(f"  - {f}")
    raise SystemExit(1)
print("\nReliable/less-reliable split + pre-PAML ORF gate tests passed.")
