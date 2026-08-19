#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FAILURES: list[str] = []


def check(cond, msg):
    if cond:
        print(f"  ok   {msg}")
    else:
        print(f"  FAIL {msg}")
        FAILURES.append(msg)


def read_tsv(path: Path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


def make_step2(work: Path, gene="T"):
    r2 = work / f"results_02/{gene}"
    r2.mkdir(parents=True, exist_ok=True)
    (r2 / f"02_{gene}.codon_for_paml.phy").write_text(
        "4 12\n"
        "A ATGAAACCCGGG\n"
        "B ATGAAACCCGGG\n"
        "C ATGAAACCCGGG\n"
        "D ATGAAACCCGGG\n"
    )
    (r2 / f"02_{gene}.tree_for_asr.nwk").write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
    return r2


def write_mock_codeml(path: Path, complete: bool, exit_code: int):
    # Four tips, but PAML itself declares THREE ancestral node records (5..7).
    # Node 5 is the degree-2 Newick serialization/root vertex; 6 and 7 are the
    # two degree-3 biological unrooted internodes.  Pensieve must validate all
    # three declared marginal records, then exclude only node 5 from biological
    # root inference downstream.
    marginal7 = "node #7 ATGAAACCCGGG\n" if complete else ""
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "from pathlib import Path\n"
        "rst = \"\"\"tree with node labels for Rod Page's TreeView\n"
        "((1_A:0.1,2_B:0.1)6:0.1,(3_C:0.1,4_D:0.1)7:0.1)5;\n"
        "Nodes 5 to 7 are ancestral\n"
        "Unreliable at sites with alignment gaps\n\n"
        "(1) Marginal reconstruction of ancestral sequences\n"
        "node #5 ATGAAACCCGGG\n"
        "node #6 ATGAAACCCGGG\n"
        f"{marginal7}"
        "Overall accuracy of the reconstruction\n"
        "(2) Joint reconstruction of ancestral sequences\n"
        "node #5 ATGAAACCCGGG\n"
        "node #6 ATGAAACCCGGG\n"
        "node #7 ATGAAACCCGGG\n"
        "\"\"\"\n"
        "Path('rst').write_text(rst)\n"
        "Path('codon_asr.out').write_text('mock codeml output\\n')\n"
        "print('Joint reconstruction.')\n"
        "print(' 50000000 bytes for conP, adjusted')\n"
        f"sys.exit({exit_code})\n"
    )
    path.chmod(0o755)


def write_slow_mock_codeml(path: Path, hang_seconds: int):
    # Simulates a real cluster failure mode: codeml finishes and flushes the
    # complete marginal ASR section to rst, prints "Joint reconstruction.",
    # then spends a long time (in reality: a large conP/space allocation that
    # can OOM-kill the job) before ever exiting. Pensieve must not sit through
    # that; it should stop codeml as soon as the marginal section is complete
    # and the Joint marker is on disk.
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import sys, time\n"
        "from pathlib import Path\n"
        "Path('rst').write_text(\n"
        "    \"tree with node labels for Rod Page's TreeView\\n\"\n"
        "    \"((1_A:0.1,2_B:0.1)6:0.1,(3_C:0.1,4_D:0.1)7:0.1)5;\\n\"\n"
        "    \"Nodes 5 to 7 are ancestral\\n\"\n"
        "    \"Unreliable at sites with alignment gaps\\n\\n\"\n"
        "    \"(1) Marginal reconstruction of ancestral sequences\\n\"\n"
        "    \"node #5 ATG AAA CCC GGG\\n\"\n"
        "    \"node #6 ATG AAA CCC GGG\\n\"\n"
        "    \"node #7 ATG AAA CCC GGG\\n\"\n"
        "    \"Overall accuracy of the reconstruction\\n\"\n"
        "    \"(2) Joint reconstruction of ancestral sequences\\n\"\n"
        ")\n"
        "Path('codon_asr.out').write_text('mock codeml output\\n')\n"
        "print('Joint reconstruction.', flush=True)\n"
        f"time.sleep({hang_seconds})\n"
        "print(' 999999999 bytes for conP, adjusted')\n"
        "sys.exit(137)\n"
    )
    path.chmod(0o755)


def test_slow_codeml_stopped_before_joint_reconstruction():
    print("codeml stopped before a slow/OOM-prone joint reconstruction")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); work = tmp / "work"; work.mkdir(); bindir = tmp / "bin"; bindir.mkdir()
        r2 = make_step2(work)
        hang_seconds = 60
        write_slow_mock_codeml(bindir / "codeml", hang_seconds=hang_seconds)
        env = os.environ.copy(); env["PATH"] = f"{bindir}:{env['PATH']}"
        env["PENSIEVE_ASR_POLL_SECONDS"] = "1"
        start = time.monotonic()
        r = subprocess.run([
            str(ROOT / "scripts/02_run_asr_backends.sh"),
            "--gene", "T", "--workdir", str(work)
        ], capture_output=True, text=True, env=env, timeout=hang_seconds)
        elapsed = time.monotonic() - start
        check(r.returncode == 0,
              f"backend succeeds instead of waiting out the joint-reconstruction hang ({r.stderr[-220:]})")
        check(elapsed < hang_seconds / 2,
              f"Pensieve stops codeml shortly after the marginal section completes rather than waiting "
              f"the full {hang_seconds}s (elapsed={elapsed:.1f}s)")
        status = read_tsv(r2 / "02_T.codeml_run_status.tsv")
        check(status and status[0]["marginal_asr_validated"] == "True"
              and status[0]["pipeline_action"] == "STOPPED_BEFORE_JOINT_RECONSTRUCTION_WITH_VALIDATED_MARGINAL_ASR",
              "run status explicitly records that codeml was stopped deliberately before joint reconstruction")
        check("stopping codeml now" in r.stderr,
              "user receives an explicit message that Pensieve stopped codeml deliberately")
        valid = read_tsv(r2 / "02_T.paml_marginal_validation.tsv")
        check(valid and valid[0]["marginal_asr_valid"] == "True"
              and valid[0]["observed_internal_sequences"] == "3",
              "codon-spaced marginal sequences (e.g. 'ATG AAA CCC GGG') are still parsed correctly")


def test_reference_free_contract():
    print("reference-free contract")
    s00 = (ROOT / "scripts/00_prune_and_check_orf.py").read_text()
    s01 = (ROOT / "scripts/01_run_macse_and_extract_events.py").read_text()
    runner = (ROOT / "scripts/run_one_gene_00_to_04.sh").read_text()
    check("--reference" not in s00, "Step 00 has no --reference option")
    check("--reference-species" not in s01, "Step 01 has no --reference-species option")
    check("reference_info.tsv" not in s00 and "reference_sequence.fasta" not in s00,
          "Step 00 creates no reference-info/reference-sequence files")
    check("relative_to_reference" not in s01,
          "Step 01 creates no reference-relative event table")
    check("primary_codon_alignment_native.fasta" in runner and "canonical_alignment.fasta" in runner,
          "canonical alignment is the alignment-coordinate file exported to final_results")

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        fa = tmp / "T.fa"
        fa.write_text("".join(f">{x}\nATGAAACCCGGG\n" for x in "ABCD"))
        tree = tmp / "tree.nwk"; tree.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
        out = tmp / "r0"
        r = subprocess.run([
            sys.executable, str(ROOT / "scripts/00_prune_and_check_orf.py"),
            "--gene", "T", "--fasta", str(fa), "--tree", str(tree), "--outdir", str(out)
        ], capture_output=True, text=True)
        check(r.returncode == 0, "Step 00 runs without any reference selection")
        names = {x.name for x in out.iterdir()}
        check(not any("reference" in x.lower() for x in names),
              "Step 00 output directory contains no reference artifacts")


def test_legacy_reference_purge():
    print("legacy reference artifact purge")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); work = tmp / "work"; work.mkdir(); bindir = tmp / "bin"; bindir.mkdir()
        for d in [
            work / "results_00/T", work / "results_01/T",
            work / "final_results/T/important_output", work / "final_results/T/supporting_files"
        ]:
            d.mkdir(parents=True, exist_ok=True)
        legacy = [
            work / "results_00/T/00_T.reference_info.tsv",
            work / "results_00/T/00_T.reference_sequence.fasta",
            work / "results_01/T/01_T.macse_indels_relative_to_reference.tsv",
            work / "final_results/T/important_output/T.reference_info.tsv",
            work / "final_results/T/important_output/T.reference_sequence.fasta",
            work / "final_results/T/supporting_files/T.macse_indels_relative_to_reference.tsv",
        ]
        for f in legacy:
            f.write_text("legacy\n")

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
    elif line.strip(): chunks.append(line.strip())
if name is not None: records.append((name,''.join(chunks)))
with outnt.open('w') as n, outaa.open('w') as a:
    for name,s in records:
        n.write(f'>{name}\n{s}\n')
        a.write(f'>{name}\n'+('X'*(len(s)//3))+'\n')
''')
        mock.chmod(0o755)
        fa = tmp / "T.fa"; fa.write_text("".join(f">{x}\nATGAAACCCGGG\n" for x in "ABCD"))
        tree = tmp / "tree.nwk"; tree.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
        env = os.environ.copy(); env["PATH"] = f"{bindir}:{env['PATH']}"
        r = subprocess.run([
            str(ROOT / "scripts/run_one_gene_00_to_04.sh"),
            "--gene", "T", "--fasta", str(fa), "--tree", str(tree),
            "--workdir", str(work), "--run_up_to", "diagnostics"
        ], capture_output=True, text=True, env=env)
        check(r.returncode == 0, f"diagnostics run used to purge legacy files ({r.stderr[-120:]})")
        check(all(not f.exists() for f in legacy),
              "known reference/reference-relative files from older workdirs are physically removed")


def test_nonzero_codeml_with_complete_marginal_continues():
    print("safe codeml non-zero handling")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); work = tmp / "work"; work.mkdir(); bindir = tmp / "bin"; bindir.mkdir()
        r2 = make_step2(work)
        write_mock_codeml(bindir / "codeml", complete=True, exit_code=1)
        env = os.environ.copy(); env["PATH"] = f"{bindir}:{env['PATH']}"
        r = subprocess.run([
            str(ROOT / "scripts/02_run_asr_backends.sh"),
            "--gene", "T", "--workdir", str(work)
        ], capture_output=True, text=True, env=env)
        check(r.returncode == 0,
              f"backend continues after codeml exit 1 when marginal ASR is complete ({r.stderr[-220:]})")
        status = read_tsv(r2 / "02_T.codeml_run_status.tsv")
        valid = read_tsv(r2 / "02_T.paml_marginal_validation.tsv")
        check(status and status[0]["codeml_exit_code"] == "1"
              and status[0]["marginal_asr_validated"] == "True"
              and status[0]["pipeline_action"] == "CONTINUE_WITH_VALIDATED_MARGINAL_ASR",
              "recovery decision explicitly records non-zero codeml + validated marginal ASR")
        check(valid and valid[0]["marginal_asr_valid"] == "True"
              and valid[0]["declared_ancestral_start"] == "5"
              and valid[0]["declared_ancestral_end"] == "7"
              and valid[0]["observed_internal_sequences"] == "3"
              and valid[0]["expected_internal_sequences"] == "3",
              "validator uses PAML's own Nodes 5 to 7 declaration and requires all three marginal sequences")
        check("continuing with the validated marginal ASR" in r.stderr,
              "user receives an explicit warning rather than a silent recovery")


def test_realistic_rst_amino_acid_section_not_treated_as_duplicate():
    # Regression for a real bug only found by running genuine codeml (not a
    # synthetic mock): real PAML rst files put "Overall accuracy of the N
    # ancestral sequences:" and an "Amino acid sequences inferred by codonml."
    # block -- itself containing "Node #N  <one-letter AA string>" lines --
    # between the marginal DNA "node #N" records and "(2) Joint
    # reconstruction". The AA lines match the same node-record pattern as the
    # real DNA records, so a parser that doesn't stop at the accuracy/AA
    # markers double-counts every node as "duplicate". This fixture mirrors
    # that real structure (trimmed) for a 3-node case (nodes 5..7).
    print("realistic rst: accuracy + amino-acid section is not mistaken for duplicate marginal records")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        phylip = tmp / "x.phy"
        phylip.write_text("4 12\nA ATGAAACCCGGG\nB ATGAAACCCGGG\nC ATGAAACCCGGG\nD ATGAAACCCGGG\n")
        rst = tmp / "rst"
        rst.write_text(
            "tree with node labels for Rod Page's TreeView\n"
            "((1_A:0.1,2_B:0.1)6:0.1,(3_C:0.1,4_D:0.1)7:0.1)5;\n"
            "Nodes 5 to 7 are ancestral\n\n"
            "(1) Marginal reconstruction of ancestral sequences\n"
            "(eqn. 4 in Yang et al. 1995 Genetics 141:1641-1650).\n\n"
            "Prob of best state at each node, listed by site\n\n"
            "   1      1   ATG (M) ATG (M) : ATG M 1.000\n\n"
            "Summary of changes along branches.\n\n"
            "Branch 1:    5..6  (n= 0.0 s= 0.0)\n\n\n"
            "List of extant and reconstructed sequences\n\n"
            "     7     12\n\n"
            "A                 ATG AAA CCC GGG \n"
            "B                 ATG AAA CCC GGG \n"
            "C                 ATG AAA CCC GGG \n"
            "D                 ATG AAA CCC GGG \n"
            "node #5           ATG AAA CCC GGG \n"
            "node #6           ATG AAA CCC GGG \n"
            "node #7           ATG AAA CCC GGG \n\n\n"
            "Overall accuracy of the 3 ancestral sequences:\n"
            "  1.00000  1.00000  1.00000\n"
            "for a site.\n\n"
            "  1.00000  1.00000  1.00000\n"
            "for the sequence.\n\n\n"
            "Amino acid sequences inferred by codonml.\n\n"
            "Node #5           MKPG\n"
            "Node #6           MKPG\n"
            "Node #7           MKPG\n\n\n"
            "Counts of changes at sites, listed by site\n\n"
            "   1 (S N:   0.000  3.000 Sd Nd:    0.0   0.0)\n\n\n"
            "(2) Joint reconstruction of ancestral sequences\n"
            "node #5           ATG AAA CCC GGG \n"
            "node #6           ATG AAA CCC GGG \n"
            "node #7           ATG AAA CCC GGG \n"
        )
        status = tmp / "status.tsv"
        r = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "02_validate_paml_marginal.py"),
             "--rst", str(rst), "--phylip", str(phylip), "--status-out", str(status)],
            capture_output=True, text=True,
        )
        check(r.returncode == 0,
              f"validator succeeds against a realistic rst with accuracy/AA sections ({r.stdout[-300:]}{r.stderr[-300:]})")
        row = read_tsv(status)[0] if status.exists() else {}
        check(row.get("marginal_asr_valid") == "True"
              and row.get("observed_internal_sequences") == "3"
              and row.get("expected_internal_sequences") == "3",
              "the accuracy/amino-acid block is not double-counted as duplicate marginal records")


def test_103_tip_rst_declared_range_regression():
    print("103-tip rst-declared ancestral range regression")
    with tempfile.TemporaryDirectory() as td:
        tmp=Path(td)
        phylip=tmp/"x.phy"
        phylip.write_text("103 3\n" + "".join(f"S{i:03d} ATG\n" for i in range(1,104)))
        rst=tmp/"rst"
        marginal="".join(f"node #{i} ATG\n" for i in range(104,206))
        joint="".join(f"node #{i} ATG\n" for i in range(104,206))
        rst.write_text(
            "Nodes 104 to 205 are ancestral\n"
            "Unreliable at sites with alignment gaps\n\n"
            "(1) Marginal reconstruction of ancestral sequences\n" + marginal +
            "(2) Joint reconstruction of ancestral sequences\n" + joint
        )
        status=tmp/"status.tsv"
        r=subprocess.run([sys.executable,str(ROOT/"scripts"/"02_validate_paml_marginal.py"),
                          "--rst",str(rst),"--phylip",str(phylip),"--status-out",str(status)],
                         capture_output=True,text=True)
        check(r.returncode==0, f"103-tip validator accepts PAML's declared 102 ancestors ({r.stdout[-160:]})")
        row=read_tsv(status)[0]
        check(row["expected_internal_sequences"]=="102" and row["observed_internal_sequences"]=="102",
              "Nodes 104..205 gives exactly 102 expected and observed marginal sequences")
        count=sum(1 for line in rst.read_text().splitlines() if "node #" in line)
        check(count==204, "fixture mirrors real rst: 102 marginal + 102 joint node records = 204 total")


def test_nonzero_codeml_with_incomplete_marginal_fails():
    print("incomplete marginal ASR remains fatal")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td); work = tmp / "work"; work.mkdir(); bindir = tmp / "bin"; bindir.mkdir()
        r2 = make_step2(work)
        write_mock_codeml(bindir / "codeml", complete=False, exit_code=1)
        env = os.environ.copy(); env["PATH"] = f"{bindir}:{env['PATH']}"
        r = subprocess.run([
            str(ROOT / "scripts/02_run_asr_backends.sh"),
            "--gene", "T", "--workdir", str(work)
        ], capture_output=True, text=True, env=env)
        check(r.returncode != 0, "backend still fails when the marginal ASR is incomplete")
        status = read_tsv(r2 / "02_T.codeml_run_status.tsv")
        check(status and status[0]["marginal_asr_validated"] == "False"
              and status[0]["pipeline_action"] == "FAIL",
              "unsafe/incomplete marginal reconstruction cannot be accepted")


def main():
    test_reference_free_contract()
    test_legacy_reference_purge()
    test_nonzero_codeml_with_complete_marginal_continues()
    test_realistic_rst_amino_acid_section_not_treated_as_duplicate()
    test_103_tip_rst_declared_range_regression()
    test_nonzero_codeml_with_incomplete_marginal_fails()
    test_slow_codeml_stopped_before_joint_reconstruction()
    if FAILURES:
        print(f"\n{len(FAILURES)} failure(s)")
        raise SystemExit(1)
    print("\nPAML-exit/reference-free tests passed.")


if __name__ == "__main__":
    main()
