#!/usr/bin/env python3
"""Step-00 input validation: every malformed input must fail loudly, and a thin
species overlap must warn without stopping.

Each of these cases previously either crashed with an opaque traceback or --
worse -- silently discarded data and produced a confident-looking result:

  * duplicate FASTA names were collapsed by a dict comprehension, so extra
    records vanished with no warning;
  * duplicate tree tips were collapsed by a set, for the same reason;
  * a protein FASTA was shredded by clean_seq() (which keeps only ACGTN-) into
    short nonsense and analysed as if it were a CDS.

Fast: step 00 only, no MACSE and no codeml.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEP00 = ROOT / "scripts" / "00_prune_and_check_orf.py"

failures = []


def check(cond, msg):
    print(f"  {'ok  ' if cond else 'FAIL'} {msg}")
    if not cond:
        failures.append(msg)


def cds(n_codons, seed=0):
    """A clean in-frame CDS: ATG + filler + no internal stop."""
    fill = ["GCT", "TTC", "AAG", "CCA", "ACC", "GGT"]
    body = "".join(fill[(seed + i) % len(fill)] for i in range(n_codons - 1))
    return "ATG" + body


def run_step00(tmp, gene, fasta_text, tree_text):
    d = Path(tmp) / gene
    d.mkdir(parents=True, exist_ok=True)
    fa, nw, out = d / "in.fas", d / "in.nwk", d / "out"
    fa.write_text(fasta_text)
    nw.write_text(tree_text)
    proc = subprocess.run(
        [sys.executable, str(STEP00), "--gene", gene, "--fasta", str(fa),
         "--tree", str(nw), "--outdir", str(out), "--alignment-mode", "perform"],
        capture_output=True, text=True)
    log = out / f"00_{gene}.input_validation.log"
    return proc, (log.read_text() if log.exists() else "")


GOOD_TREE = "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);"
GOOD_FASTA = "".join(f">{n}\n{cds(20, i)}\n" for i, n in enumerate("ABCD"))


def main():
    with tempfile.TemporaryDirectory() as tmp:
        print("step 00 input validation")

        # --- control: a well-formed input must succeed -----------------------
        p, log = run_step00(tmp, "GOOD", GOOD_FASTA, GOOD_TREE)
        check(p.returncode == 0, "well-formed input succeeds")
        check("[OK] 4 species shared" in log, "control logs the shared-species count")

        # --- halts -----------------------------------------------------------
        dup_fa = GOOD_FASTA + f">A\n{cds(20, 9)}\n"
        p, log = run_step00(tmp, "DUPFASTA", dup_fa, GOOD_TREE)
        check(p.returncode != 0, "duplicate FASTA names halt")
        check("duplicate sequence names" in log and "A (x2)" in log,
              "duplicate FASTA names are named in the log")

        dup_tree = "(((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2):0.1,A:0.4);"
        p, log = run_step00(tmp, "DUPTREE", GOOD_FASTA, dup_tree)
        check(p.returncode != 0, "duplicate tree tips halt")
        check("appear more than once" in log and "A (x2)" in log,
              "duplicate tree tips are named in the log")

        p, log = run_step00(tmp, "EMPTYTREE", GOOD_FASTA, "   \n\n")
        check(p.returncode != 0, "empty tree halts")
        check("--tree is empty" in log, "empty tree says so explicitly")

        p, log = run_step00(tmp, "BADTREE", GOOD_FASTA, "#NEXUS\nTREE t = ((A,B),C\n")
        check(p.returncode != 0, "unparseable tree halts")
        check("could not be parsed as a Newick tree" in log, "bad tree names the parse failure")

        prot = "".join(f">{n}\nMKPQEFILWYVGAT\n" for n in "ABCD")
        p, log = run_step00(tmp, "PROTEIN", prot, GOOD_TREE)
        check(p.returncode != 0, "amino-acid FASTA halts")
        check("does not look like NUCLEOTIDE" in log, "protein FASTA says nucleotide is required")

        renamed = "".join(f">X_{n}\n{cds(20, i)}\n" for i, n in enumerate("ABCD"))
        p, log = run_step00(tmp, "NOMATCH", renamed, GOOD_TREE)
        check(p.returncode != 0, "zero species overlap halts")
        check("no species name is shared" in log, "zero overlap explains the name mismatch")

        # --- warn and continue ------------------------------------------------
        # 4 tree tips, 4 FASTA records, only 1 shared -> 25%
        thin = f">A\n{cds(20)}\n" + "".join(f">Z_{n}\n{cds(20, i)}\n" for i, n in enumerate("BCD"))
        p, log = run_step00(tmp, "THIN", thin, GOOD_TREE)
        check(p.returncode == 0, "thin species overlap continues rather than halting")
        check("[WARNING]" in log and "below 50%" in log, "thin overlap warns in the log")

        # --- all-complete / all-incomplete notes -------------------------------
        p, log = run_step00(tmp, "ALLOK", GOOD_FASTA, GOOD_TREE)
        check("every sequence is a complete ORF" in log, "all-complete case is noted")

        broken = "".join(f">{n}\n{'ATGTAACCC' + cds(10, i)}\n" for i, n in enumerate("ABCD"))
        p, log = run_step00(tmp, "ALLBROKEN", broken, GOOD_TREE)
        check(p.returncode == 0, "all-incomplete input still succeeds")
        check("no sequence has a complete ORF" in log, "all-incomplete case is noted")

    print()
    if failures:
        print(f"{len(failures)} failure(s)")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print("Input validation test passed.")


if __name__ == "__main__":
    main()
