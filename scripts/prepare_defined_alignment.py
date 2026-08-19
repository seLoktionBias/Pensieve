#!/usr/bin/env python3
"""Opt-in preprocessing to make a curated alignment ready for --alignment defined.

Pensieve's core stays deliberately strict: 02_prepare_asr_inputs.py refuses any
defined alignment whose column count is not a multiple of three and never adds,
removes, or reorders columns on its own (see INFERENCE_SPEC / CHANGELOG). That is
correct for the authoritative-alignment contract, but some curated alignments
arrive with two purely mechanical issues that are safe to resolve WITHOUT any
biological judgement, and this standalone helper resolves exactly those two --
and nothing else -- writing a new FASTA the user can then pass to --alignment
defined. It never edits real residues and never touches internal columns.

1. Assembly-gap masking (optional, --mask-regions):
   Per-species regions listed as `species<TAB>start<TAB>end` (1-based, inclusive,
   in the INPUT alignment's own columns) mark stretches the collaborator flags as
   unreliable assembly gaps. Within each region, only gap characters ('-') are
   rewritten to 'N' (explicit "unknown"), so codeml (cleandata=0) treats them as
   ambiguous rather than as confident deletions. Real bases are never changed.
   This reproduces the already_aligned_bat_genes/replace_alignment_gaps.py step.

2. Trailing/leading all-gap padding trim (frame fix):
   A column that is a gap in EVERY sequence carries zero information and, when it
   sits at the very start or end of the alignment (before/after all real
   content), can be removed without altering any residue's homology or shifting
   any codon boundary of the real content. INTERNAL all-gap columns are left
   untouched: they are legitimate codon-frame insertion columns, and removing
   them WOULD misframe everything downstream (verified on real CNGB3 data: after
   trimming only the 8 trailing all-gap columns the alignment reads as a clean
   codon frame for the intact-ORF species; removing the internal all-gap columns
   instead introduced spurious in-frame stops in every intact species). Only pure
   '-' columns are trimmed, so masking-introduced 'N' is never discarded.

After these two steps the alignment length must be a multiple of three. If it is
not, the residual frame problem is genuine (real, non-gap content is out of
codon frame) and CANNOT be fixed mechanically without altering data -- the script
refuses and reports the residual, rather than silently trimming real residues.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO


def read_alignment(path: Path):
    records, order = {}, []
    for record in SeqIO.parse(str(path), "fasta"):
        name = record.id.split()[0]
        if name in records:
            raise SystemExit(f"{path}: duplicate FASTA identifier {name}")
        records[name] = list(str(record.seq).upper().replace("U", "T"))
        order.append(name)
    if not records:
        raise SystemExit(f"{path}: no FASTA records")
    lengths = {len(records[name]) for name in order}
    if len(lengths) != 1:
        raise SystemExit(f"{path}: sequences are not equal length: {sorted(lengths)}")
    return records, order, lengths.pop()


def load_mask_regions(path: Path):
    regions: dict[str, list[tuple[int, int]]] = {}
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 3:
                raise SystemExit(f"{path}: expected 'species start end' per line, got: {line!r}")
            species, start, end = parts[0], int(parts[1]), int(parts[2])
            if start > end:
                start, end = end, start
            regions.setdefault(species, []).append((start, end))
    return regions


def apply_masking(records, order, length, regions):
    """Rewrite gaps ('-') to 'N' inside each species' flagged regions only."""
    masked = 0
    for species, spans in regions.items():
        if species not in records:
            print(f"[warn] mask-regions names species not in alignment, skipped: {species}")
            continue
        seq = records[species]
        for start, end in spans:
            if start < 1 or end > length:
                print(f"[warn] {species} region {start}-{end} outside alignment length {length}; clamped")
            for i in range(max(1, start) - 1, min(length, end)):
                if seq[i] == "-":
                    seq[i] = "N"
                    masked += 1
    return masked


def trim_terminal_allgap(records, order, length):
    """Remove leading and trailing columns that are '-' in every sequence."""
    def col_all_gap(i):
        return all(records[name][i] == "-" for name in order)

    lead = 0
    while lead < length and col_all_gap(lead):
        lead += 1
    trail = 0
    while trail < length - lead and col_all_gap(length - 1 - trail):
        trail += 1
    if lead or trail:
        keep = range(lead, length - trail)
        for name in order:
            records[name] = records[name][lead:length - trail]
    return lead, trail, length - lead - trail


def write_alignment(records, order, path: Path, width=60):
    with path.open("w") as handle:
        for name in order:
            seq = "".join(records[name])
            handle.write(f">{name}\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i:i + width] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="infile", required=True, help="input codon alignment FASTA")
    ap.add_argument("--out", dest="outfile", required=True, help="preprocessed FASTA for --alignment defined")
    ap.add_argument("--mask-regions", default=None,
                    help="optional TSV of species<TAB>start<TAB>end assembly-gap regions (1-based, inclusive)")
    args = ap.parse_args()

    records, order, length = read_alignment(Path(args.infile))
    print(f"[info] input: {len(order)} sequences x {length} columns (length mod 3 = {length % 3})")

    if args.mask_regions:
        regions = load_mask_regions(Path(args.mask_regions))
        n_masked = apply_masking(records, order, length, regions)
        print(f"[info] assembly-gap masking: {n_masked} gap position(s) rewritten to N "
              f"across {len(regions)} species")

    lead, trail, new_length = trim_terminal_allgap(records, order, length)
    print(f"[info] trimmed all-gap padding: {lead} leading + {trail} trailing column(s) "
          f"-> {new_length} columns (length mod 3 = {new_length % 3})")

    if new_length % 3 != 0:
        raise SystemExit(
            f"[error] after masking and all-gap padding trim the alignment is {new_length} columns "
            f"(mod 3 = {new_length % 3}), still not a codon-frame length. The residual offset is real, "
            "non-gap content out of frame and cannot be fixed mechanically without altering sequence "
            "data. Reframe the alignment manually (codon-aware) before using --alignment defined."
        )

    write_alignment(records, order, Path(args.outfile))
    print(f"[done] wrote {args.outfile}: {len(order)} sequences x {new_length} columns, ready for --alignment defined")


if __name__ == "__main__":
    main()
