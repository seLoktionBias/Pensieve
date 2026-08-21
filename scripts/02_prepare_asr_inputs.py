#!/usr/bin/env python3
"""Prepare one authoritative alignment and synchronized downstream views.

Pensieve v3.30 design
---------------------
* --alignment perform: the MACSE NT alignment is the canonical alignment.
  MACSE is the only alignment engine; no second aligner or coordinate
  projection is performed.
* --alignment defined: the pruned user alignment is canonical and its columns
  are never inserted, deleted, or reordered by Pensieve.
* The native view preserves observed residues/gaps and renders MACSE '!' as '-'
  (removing MACSE's frame-restoration placeholder).  This does NOT mean each '!'
  is itself a deletion; insertion/deletion direction is inferred later from the
  complete aligned occupancy pattern and the phylogeny.
* The PAML-safe view has exactly the same columns. MACSE '!' is rendered as 'N'
  and exact stop codons are masked to NNN.  cleandata=0 in codeml retains these
  ambiguous codons.
* Raw premature STOPs are mapped to canonical alignment coordinates before they
  are masked, so shared nonsense events cannot disappear during preparation.
"""
from __future__ import annotations

import argparse
import csv
import difflib
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import Phylo, SeqIO
from Bio.Data import CodonTable

DNA = set("ACGT")
DNA_N = set("ACGTN")
ALLOWED_DEFINED = set("ACGTN-")
ALLOWED_MACSE = set("ACGTN-!")
STOPS = {"TAA", "TAG", "TGA"}


def read_tsv(path: Path) -> List[dict]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(rows: List[dict], path: Path, header: List[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if header is None:
        header = []
        for row in rows:
            for key in row:
                if key not in header:
                    header.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=header, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> Tuple[Dict[str, str], List[str]]:
    records, order = {}, []
    for record in SeqIO.parse(str(path), "fasta"):
        name = record.id.split()[0]
        if name in records:
            raise SystemExit(f"{path}: duplicate FASTA identifier {name}")
        seq = str(record.seq).upper().replace("U", "T")
        records[name] = seq
        order.append(name)
    if not records:
        raise SystemExit(f"{path}: no FASTA records")
    return records, order


def write_fasta(records: Dict[str, str], path: Path, order: List[str]) -> None:
    with path.open("w") as handle:
        for name in order:
            seq = records[name]
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def write_phylip(records: Dict[str, str], path: Path, order: List[str]) -> None:
    if not order:
        raise SystemExit("Cannot write PHYLIP with no sequences")
    lengths = {len(records[name]) for name in order}
    if len(lengths) != 1:
        raise SystemExit(f"Cannot write PHYLIP: unequal lengths {sorted(lengths)}")
    length = lengths.pop()
    with path.open("w") as handle:
        handle.write(f" {len(order)} {length}\n")
        for name in order:
            handle.write(f"{name}  {records[name]}\n")


def validate_equal_alignment(records: Dict[str, str], order: List[str], allowed: set[str], label: str) -> int:
    if not order:
        raise SystemExit(f"{label}: no sequences")
    lengths = {len(records[name]) for name in order}
    if len(lengths) != 1:
        raise SystemExit(f"{label}: sequences have unequal lengths: {sorted(lengths)}")
    length = lengths.pop()
    bad = []
    for name in order:
        chars = sorted(set(records[name]) - allowed)
        if chars:
            bad.append(f"{name}:{''.join(chars)}")
    if bad:
        raise SystemExit(f"{label}: unsupported sequence characters: {'; '.join(bad[:20])}")
    if length % 3 != 0:
        raise SystemExit(
            f"{label}: alignment length {length} is not divisible by 3. Pensieve v3.30 will not "
            "silently add/remove alignment columns to manufacture a codon boundary."
        )
    return length


def build_raw_to_alignment_map(aligned: str, raw_seq: str) -> Dict[int, int]:
    """Build a full 1-based-raw-position -> 1-based-alignment-column map for
    one species' own row, once, for repeated lookup.

    MACSE never deletes or reorders a real input character, but it can
    insert synthetic characters not present in the raw input at all -- not
    only its own documented '!' partial-codon placeholder, but (confirmed on
    real data: GUCA1C/Desmodus_rotundus) occasionally a literal ambiguity
    code like 'N' standing in for a position it could not resolve, with no
    '!' involved, and sometimes more than once in the same row. Skipping
    only '-' and '!' therefore silently miscounts whenever MACSE emits one
    of these. A naive single-pass greedy character match (advance on match,
    otherwise skip) also breaks down as soon as there is more than one such
    insertion in the same row -- it has no way to recover once a real
    character fails to match by coincidence, and can stall with no further
    matches for the rest of the sequence (confirmed on real data: this
    exact failure mode blocked GUCA1C/Desmodus_rotundus's STOP-coordinate
    mapping entirely). See build_raw_to_alignment_map() in
    01_run_macse_and_extract_events.py for the full real-data account.

    A proper sequence alignment (difflib's longest-matching-block algorithm)
    between the species' own raw sequence and the non-gap characters of its
    aligned row does not have this failure mode: every 'equal' opcode block
    is, by construction, a genuine, correctly-ordered run of real raw
    characters, however many separate synthetic insertions sit between them.
    """
    non_gap_positions = [i for i, c in enumerate(aligned) if c != "-"]
    non_gap_chars = "".join(aligned[i] for i in non_gap_positions)
    sm = difflib.SequenceMatcher(None, raw_seq.upper(), non_gap_chars.upper(), autojunk=False)
    mapping: Dict[int, int] = {}
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag != "equal":
            continue
        for k in range(i2 - i1):
            mapping[i1 + k + 1] = non_gap_positions[j1 + k] + 1
    return mapping


def safe_int(value) -> int | None:
    try:
        return int(float(str(value)))
    except Exception:
        return None


def trueish(value) -> bool:
    return str(value).strip().lower() in {"true", "t", "yes", "1"}


def codon_to_aa(codon: str) -> str:
    codon = codon.upper().replace("U", "T")
    if len(codon) != 3:
        return "X"
    if codon == "---":
        return "-"
    if any(ch not in DNA for ch in codon):
        return "X"
    if codon in STOPS:
        return "*"
    return CodonTable.unambiguous_dna_by_name["Standard"].forward_table.get(codon, "X")


def aa_alignment(nt: Dict[str, str], order: List[str]) -> Dict[str, str]:
    out = {}
    for name in order:
        seq = nt[name]
        out[name] = "".join(codon_to_aa(seq[i:i + 3]) for i in range(0, len(seq), 3))
    return out


def load_raw_stop_diagnostics(r0: Path, r1: Path, gene: str) -> List[dict]:
    """Join raw STOP calls to MACSE frame-phase diagnostics.

    v3.30 no longer classifies a STOP as a frameshift consequence merely because
    any ``!`` marker occurs upstream.  Step 01 reports the cumulative MACSE frame
    correction modulo three at the STOP itself.  A compensated frame (mod3 == 0)
    leaves the STOP eligible as an independent nonsense event.
    """
    phase = {}
    for row in read_tsv(r1 / f"01_{gene}.macse_premature_stop_masking.tsv"):
        key = (row.get("species", ""), safe_int(row.get("raw_nt_start")), safe_int(row.get("raw_nt_end")))
        shifted = row.get("frame_shifted_at_stop")
        if shifted in {None, "", "NA"}:
            # Backward-compatible fallback for an older Step-01 table.  This is
            # less informative than v3.30 but avoids misreading the field.
            shifted = row.get("masked_by_upstream_macse_frameshift_marker", False)
        phase[key] = {
            "frame_shifted_at_stop": trueish(shifted),
            "upstream_macse_marker_count": row.get("upstream_macse_marker_count", "NA"),
            "upstream_macse_frame_correction_mod3": row.get("upstream_macse_frame_correction_mod3", "NA"),
            "stop_phase_interpretation": row.get("interpretation", "NA"),
        }

    rows = []
    for row in read_tsv(r0 / f"00_{gene}.orf_failures.tsv"):
        if row.get("failure_type") != "premature_in_frame_stop":
            continue
        start = safe_int(row.get("nt_start"))
        end = safe_int(row.get("nt_end"))
        if start is None or end is None:
            continue
        key = (row.get("species", ""), start, end)
        diag = phase.get(key, {})
        rows.append({
            "gene": gene,
            "species": row.get("species", ""),
            "raw_codon_position": row.get("codon_position", "NA"),
            "raw_nt_start": start,
            "raw_nt_end": end,
            "stop_codon": str(row.get("codon", "NA")).upper(),
            "frame_shifted_at_stop": diag.get("frame_shifted_at_stop", False),
            "upstream_macse_marker_count": diag.get("upstream_macse_marker_count", "NA"),
            "upstream_macse_frame_correction_mod3": diag.get("upstream_macse_frame_correction_mod3", "NA"),
            "stop_phase_interpretation": diag.get("stop_phase_interpretation", "NA"),
        })
    return rows


def build_stop_registry(gene: str, raw_stops: List[dict], canonical_source: Dict[str, str],
                        raw_by_species: Dict[str, str], coordinate_system: str) -> List[dict]:
    """Map each raw premature STOP into canonical alignment coordinates and
    decide whether it may found an independent nonsense-mutation event.

    Four independent gates must all pass, matching the ChatGPT-review
    finding that a raw STOP alone is not enough evidence (this is the fix
    for the CNGA3/Phyllostomus discolor overfitting regression -- a single
    early frameshift there produced many downstream raw STOP triplets that
    disappear once that one frameshift is corrected, and those must not be
    reported as independent nonsense substitutions):

    1. ``independent`` -- step01's frame-phase bookkeeping (cumulative
       upstream MACSE frame correction mod 3) says the reading frame is not
       shifted at this exact STOP.
    2. ``contiguous`` -- the raw STOP's three nucleotides map to three
       CONSECUTIVE canonical alignment columns for this species (no MACSE
       '!' placeholder or gap wedged between them). A non-contiguous mapping
       means the "codon" spans a structurally disturbed region and cannot
       be trusted as a clean lesion.
    3. ``codon_frame_aligned`` -- the mapped span actually starts on a real
       codon boundary of the canonical alignment's ONE shared codon frame
       (column 1 = the first base of codon 1, column 4 = the first base of
       codon 2, and so on for every species alike -- this is the same frame
       every other codon-level operation in the pipeline, e.g. codon_to_aa/
       aa_alignment(), already assumes). Real bug, found on real CNGA3 data
       (reported directly by the user): gate 1 alone is not sufficient. A
       species can have upstream MACSE markers whose LENGTHS sum to a
       multiple of 3 (so gate 1 calls the frame "restored") while still
       shifting which RAW nucleotides fall into which codon NUMBER from
       that point on, because the raw premature-stop scanner
       (00_prune_and_check_orf.py) numbers codons by counting from raw
       position 1 with no knowledge of any correction. "Raw codon 18" and
       "the alignment's actual codon 18" can therefore be two different,
       non-overlapping spans once an odd number of markers happens to
       total a multiple of 3 in aggregate. Confirmed directly on real data:
       CNGA3/Phyllostomus_discolor's raw-scanner-reported codon 18 mapped
       to alignment columns 53-55 -- NOT a codon boundary (52 is not a
       multiple of 3) -- and reads as the last two bases of the alignment's
       real codon 18 (translates to Valine) plus the first base of its real
       codon 19 (Arginine): "TAA" purely by coincidence of an off-frame
       window, not a real stop codon. Phyllostomus_discolor's own protein,
       translated in its real, frame-aligned codons, matches its relatives
       cleanly with no disruption at this position at all.
    4. ``codon_confirmed`` -- the actual homologous codon read directly from
       the canonical alignment at those columns is still exactly
       TAA/TAG/TGA. This re-derives the codon from the FINAL alignment
       (after conserved-block/frame corrections) rather than trusting the
       RAW pre-alignment codon text unconditionally.
    """
    mapped = []
    raw_to_aln_by_species: Dict[str, Dict[int, int]] = {}
    for raw in raw_stops:
        species = raw["species"]
        seq = canonical_source.get(species)
        raw_seq = raw_by_species.get(species)
        if seq is None or raw_seq is None:
            continue
        if species not in raw_to_aln_by_species:
            raw_to_aln_by_species[species] = build_raw_to_alignment_map(seq, raw_seq)
        raw_to_aln = raw_to_aln_by_species[species]
        cols = [raw_to_aln.get(pos)
                for pos in range(raw["raw_nt_start"], raw["raw_nt_end"] + 1)]
        valid = all(col is not None for col in cols)
        contiguous = valid and cols == list(range(cols[0], cols[0] + len(cols)))
        aligned_start = min(cols) if valid else None
        aligned_end = max(cols) if valid else None
        clean_span = bool(contiguous and valid and aligned_end - aligned_start + 1 == 3)
        frame_aligned = bool(clean_span and (aligned_start - 1) % 3 == 0)
        corrected_codon = seq[aligned_start - 1:aligned_end].upper() if clean_span else None
        codon_confirmed = bool(frame_aligned and corrected_codon in STOPS)
        independent = not raw["frame_shifted_at_stop"]
        candidate = bool(independent and valid and codon_confirmed)
        if not independent:
            reason = "raw_premature_stop_in_shifted_frame_recorded_as_frameshift_consequence"
        elif not valid:
            reason = "raw_premature_stop_could_not_be_mapped_to_canonical_alignment"
        elif not clean_span:
            reason = "raw_premature_stop_non_contiguous_or_non_3bp_mapping_uncertain"
        elif not frame_aligned:
            reason = "raw_premature_stop_mapped_span_does_not_start_on_a_real_codon_boundary"
        elif not codon_confirmed:
            reason = "raw_premature_stop_corrected_homologous_codon_is_not_a_stop"
        else:
            reason = "raw_premature_stop_mapped_to_canonical_alignment_as_independent_candidate"
        mapped.append({
            "gene": gene,
            "species": species,
            "codon_position": raw["raw_codon_position"],
            "nt_start": raw["raw_nt_start"],
            "nt_end": raw["raw_nt_end"],
            "primary_alignment_start": aligned_start if aligned_start is not None else "NA",
            "primary_alignment_end": aligned_end if aligned_end is not None else "NA",
            "stop_codon": raw["stop_codon"],
            "corrected_homologous_codon": corrected_codon if corrected_codon is not None else "NA",
            "coordinate_system": coordinate_system,
            "terminal_codon": False,
            "pseudogenizing_event_candidate": candidate,
            "independent_stop_candidate": candidate,
            "mapped_columns_contiguous": contiguous,
            "codon_frame_aligned": frame_aligned,
            "codon_confirmed_stop": codon_confirmed,
            "frame_shifted_at_stop": raw["frame_shifted_at_stop"],
            "upstream_macse_marker_count": raw["upstream_macse_marker_count"],
            "upstream_macse_frame_correction_mod3": raw["upstream_macse_frame_correction_mod3"],
            "stop_phase_interpretation": raw["stop_phase_interpretation"],
            "reason": reason,
            "stop_event_key": (
                f"{aligned_start}-{aligned_end}:{raw['stop_codon']}"
                if valid else "NA"
            ),
        })
    return mapped


def scan_defined_stops(gene: str, native: Dict[str, str], order: List[str],
                       coordinate_system: str) -> List[dict]:
    """--alignment defined: detect premature stops by reading the authoritative
    codon alignment directly, with no MACSE and no raw->alignment remapping.

    The alignment length is a multiple of three (already validated), so every
    codon is simply a consecutive column triplet -- the same thing as running
    ``fold -w 3`` on each row, which never masks a TAA/TAG/TGA that straddles a
    codon boundary because there is no boundary to straddle here: column 1 is
    the first base of codon 1 for every species alike. For each species we walk
    those codons, skip its own terminal codon (the terminal stop was already
    removed by gapping in step 00, but we guard against any that remain), and
    record every in-frame TAA/TAG/TGA at its exact alignment columns as an
    independent premature-stop candidate. Because these coordinates are read
    straight off the canonical alignment, the four MACSE-era gates
    (frame-phase, contiguity, codon-frame, codon-confirmation) are all
    satisfied by construction, so each is reported as already met.
    """
    length = len(next(iter(native.values())))
    rows: List[dict] = []
    for species in order:
        seq = native[species].upper()
        last_nongap = max((i for i, ch in enumerate(seq) if ch != "-"), default=-1)
        terminal_codon_idx = last_nongap // 3 if last_nongap >= 0 else -1
        for col0 in range(0, length, 3):
            codon = seq[col0:col0 + 3]
            if codon not in STOPS:
                continue
            codon_idx = col0 // 3
            is_terminal = codon_idx == terminal_codon_idx
            start, end = col0 + 1, col0 + 3
            rows.append({
                "gene": gene,
                "species": species,
                "codon_position": codon_idx + 1,
                "nt_start": "NA",
                "nt_end": "NA",
                "primary_alignment_start": start,
                "primary_alignment_end": end,
                "stop_codon": codon,
                "corrected_homologous_codon": codon,
                "coordinate_system": coordinate_system,
                "terminal_codon": is_terminal,
                "pseudogenizing_event_candidate": (not is_terminal),
                "independent_stop_candidate": (not is_terminal),
                "mapped_columns_contiguous": True,
                "codon_frame_aligned": True,
                "codon_confirmed_stop": True,
                "frame_shifted_at_stop": "NA",
                "upstream_macse_marker_count": "NA",
                "upstream_macse_frame_correction_mod3": "NA",
                "stop_phase_interpretation": "defined_alignment_direct_codon_scan",
                "reason": (
                    "defined_alignment_terminal_stop_masked_for_paml_only" if is_terminal
                    else "defined_alignment_inframe_stop_masked_as_independent_candidate"
                ),
                "stop_event_key": f"{start}-{end}:{codon}",
            })
    return rows


def mask_paml_stops(paml: Dict[str, str], order: List[str], registry: List[dict],
                    gene: str, coordinate_system: str) -> Tuple[Dict[str, str], List[dict]]:
    """Mask exact aligned stop codons but retain a complete audit trail."""
    out = {name: list(paml[name]) for name in order}
    known_raw = {(r["species"], safe_int(r["primary_alignment_start"]), safe_int(r["primary_alignment_end"])): r
                 for r in registry}
    technical = []
    length = len(next(iter(paml.values())))
    for species in order:
        for start0 in range(0, length, 3):
            codon = "".join(out[species][start0:start0 + 3]).upper()
            if codon not in STOPS:
                continue
            start, end = start0 + 1, start0 + 3
            raw = known_raw.get((species, start, end))
            is_terminal = start0 == length - 3
            if raw is None:
                technical.append({
                    "gene": gene, "species": species, "codon_position": start0 // 3 + 1,
                    "nt_start": "NA", "nt_end": "NA",
                    "primary_alignment_start": start, "primary_alignment_end": end,
                    "stop_codon": codon, "corrected_homologous_codon": codon,
                    "coordinate_system": coordinate_system,
                    "terminal_codon": is_terminal,
                    "pseudogenizing_event_candidate": False,
                    "independent_stop_candidate": False,
                    "mapped_columns_contiguous": True,
                    "codon_frame_aligned": True,
                    "codon_confirmed_stop": True,
                    "frame_shifted_at_stop": "NA",
                    "upstream_macse_marker_count": "NA",
                    "upstream_macse_frame_correction_mod3": "NA",
                    "stop_phase_interpretation": "NA",
                    "reason": "terminal_or_unregistered_stop_masked_for_paml_only",
                    "stop_event_key": f"{start}-{end}:{codon}",
                })
            out[species][start0:start0 + 3] = list("NNN")

    # Also mask mapped raw stops when the aligned triplet is not an exact global
    # codon window. This keeps codeml from seeing a known stop while preserving
    # its original event coordinates in the registry.
    for row in registry:
        start = safe_int(row["primary_alignment_start"])
        end = safe_int(row["primary_alignment_end"])
        species = row["species"]
        if start is None or end is None:
            continue
        for i in range(start - 1, min(end, length)):
            if out[species][i] != "-":
                out[species][i] = "N"
    return {name: "".join(out[name]) for name in order}, technical


COMPLETE_ORF_VALIDATION_HEADER = [
    "gene", "species", "complete_orf", "input_length", "aligned_length", "degapped_length",
    "degapped_matches_input", "n_frameshift_markers", "n_incomplete_codon_cells",
    "first_incomplete_codon_column", "incomplete_codon_cells", "verdict", "failure_reason",
]


def validate_complete_orf_alignment(gene, canonical_source, native, order, complete_species,
                                    input_seqs, out, coordinate_system, on_violation="warn"):
    """Hard gate, run before anything is handed to PAML: a complete ORF must
    come out of alignment still intact.

    A sequence that entered step 00 with a complete ORF (ATG start, length a
    multiple of three, no internal in-frame STOP) describes a real, unbroken
    reading frame. Alignment is only allowed to insert gaps around it -- never
    to edit it, never to break its codon phase. Every complete-ORF row must
    therefore satisfy all three of:

      1. `degapped_matches_input` -- stripping the gap characters from its
         aligned row reproduces its input sequence EXACTLY. Any mismatch means
         alignment invented, dropped or altered a real nucleotide.
      2. `n_frameshift_markers == 0` -- no MACSE `!` partial-codon placeholder
         anywhere in the row. A `!` in a complete ORF is MACSE asserting a
         frameshift in a sequence that by construction has none; step 01 sets
         `-fs`/`-fs_term` prohibitively high for exactly this reason, and this
         gate confirms it worked.
      3. `n_incomplete_codon_cells == 0` -- every codon cell (columns 3k+1..3k+3)
         is either three real residues or exactly `---`. A mixed cell such as
         `AT-` means the alignment split a codon across a gap boundary, so the
         codon frame that PAML, the STOP scan and the ORF walk all assume is
         not actually there.

    `on_violation` decides what a violation costs:

      * `"warn"` (default) -- report it loudly, flag the species, and carry on.
        One bad row in one gene must not abandon a whole batch; the affected
        species are marked with a trailing `*` on both figures so the reader can
        see at a glance which reading frames were not verifiable.
      * `"stop"` -- abort before PAML. Use when a corrupted reading frame must
        never reach the ancestral reconstruction unnoticed.

    Either way `02_<gene>.complete_orf_alignment_validation.tsv` is written, on a
    clean run too, so the result is auditable rather than merely silent.
    """
    rows = []
    failures = []
    complete = set(complete_species)
    for species in order:
        aligned = native[species]
        source = canonical_source.get(species, aligned)
        degapped = aligned.replace("-", "")
        expected = input_seqs.get(species)
        n_markers = source.count("!")
        bad_cells = []
        for start in range(0, len(aligned), 3):
            cell = aligned[start:start + 3]
            if len(cell) != 3:
                bad_cells.append(start + 1)
                continue
            gaps = cell.count("-")
            if gaps not in (0, 3):
                bad_cells.append(start + 1)
        matches = expected is not None and degapped == expected
        is_complete = species in complete
        reasons = []
        if is_complete:
            if not matches:
                reasons.append("degapped_row_does_not_reproduce_input_sequence"
                               if expected is not None else "input_sequence_unavailable")
            if n_markers:
                reasons.append(f"{n_markers}_frameshift_marker(s)_in_a_complete_ORF")
            if bad_cells:
                reasons.append(f"{len(bad_cells)}_incomplete_codon_cell(s)")
        verdict = "not_applicable" if not is_complete else ("pass" if not reasons else "FAIL")
        row = {
            "gene": gene, "species": species, "complete_orf": is_complete,
            "input_length": len(expected) if expected is not None else "NA",
            "aligned_length": len(aligned), "degapped_length": len(degapped),
            "degapped_matches_input": matches,
            "n_frameshift_markers": n_markers,
            "n_incomplete_codon_cells": len(bad_cells),
            "first_incomplete_codon_column": bad_cells[0] if bad_cells else "NA",
            "incomplete_codon_cells": ",".join(str(c) for c in bad_cells[:50]) or "NA",
            "verdict": verdict,
            "failure_reason": ";".join(reasons) or "NA",
        }
        rows.append(row)
        if verdict == "FAIL":
            failures.append(row)

    report = out / f"02_{gene}.complete_orf_alignment_validation.tsv"
    write_tsv(rows, report, COMPLETE_ORF_VALIDATION_HEADER)
    n_complete = sum(1 for r in rows if r["complete_orf"])
    if failures:
        detail = "\n".join(
            f"    {r['species']}: {r['failure_reason']}" for r in failures[:20])
        more = f"\n    ... and {len(failures) - 20} more" if len(failures) > 20 else ""
        headline = (f"{gene}: {len(failures)} of {n_complete} complete-ORF sequence(s) did not survive "
                    f"alignment intact")
        if on_violation == "stop":
            raise SystemExit(
                f"[ERROR] {headline}. Pensieve stops here rather than handing a corrupted reading "
                f"frame to PAML (--on-complete-orf-violation stop).\n{detail}{more}\n"
                f"    Diagnostic report: {report}\n"
                f"    Coordinate system: {coordinate_system}"
            )
        print(f"[WARN] {headline}. Continuing (--on-complete-orf-violation warn); these species are "
              f"marked with a trailing '*' on both figures and their reading frame downstream is "
              f"NOT verified.\n{detail}{more}\n"
              f"    Diagnostic report: {report}", file=sys.stderr)
    else:
        print(f"Complete-ORF alignment validation for {gene}: {n_complete}/{n_complete} complete ORFs "
              f"reproduce their input exactly, carry no frameshift markers, and sit on whole codon cells "
              f"({report.name})")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Pensieve's canonical native and PAML-safe alignments.")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--results01-dir", required=True)
    parser.add_argument("--results00-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--alignment-mode", choices=["perform", "defined"], default="perform")
    parser.add_argument("--on-complete-orf-violation", choices=["stop", "warn"], default="warn",
                        help="What to do when a complete-ORF sequence does not survive alignment intact "
                             "(see validate_complete_orf_alignment). 'warn' (default) reports it, flags "
                             "the species with a trailing '*' on both figures, and continues -- one bad "
                             "row must not abandon a whole batch. 'stop' aborts before PAML. The "
                             "diagnostic report 02_<gene>.complete_orf_alignment_validation.tsv is "
                             "written either way.")
    args = parser.parse_args()

    gene = args.gene
    r0 = Path(args.results00_dir).resolve()
    r1 = Path(args.results01_dir).resolve()
    out = Path(args.outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    user, user_order = read_fasta(r0 / f"00_{gene}.common_species.fasta")
    tree = Phylo.read(str(r0 / f"00_{gene}.common_species.tree"), "newick")
    def terminals_top_to_bottom(clade):
        if not clade.clades:
            return [clade]
        tips = []
        for child in reversed(clade.clades):
            tips.extend(terminals_top_to_bottom(child))
        return tips
    tree_order = [tip.name for tip in terminals_top_to_bottom(tree.root)]
    order = [name for name in tree_order if name in user]
    if set(order) != set(user):
        raise SystemExit("Step02 species mismatch between pruned FASTA and pruned tree")

    if args.alignment_mode == "perform":
        canonical_nt_name = f"01_{gene}.macse_NT.fasta"
        macse, macse_order = read_fasta(r1 / canonical_nt_name)
        if set(macse) != set(user):
            raise SystemExit(
                f"MACSE output species differ from pruned input: missing={sorted(set(user)-set(macse))[:10]}, "
                f"extra={sorted(set(macse)-set(user))[:10]}"
            )
        validate_equal_alignment(macse, order, ALLOWED_MACSE, "MACSE alignment")
        canonical_source = {name: macse[name] for name in order}
        native = {name: macse[name].replace("!", "-") for name in order}
        paml_pre = {name: macse[name].replace("!", "N") for name in order}
        coordinate_system = "macse_codon_alignment"
        source_path = r1 / canonical_nt_name
    else:
        validate_equal_alignment(user, order, ALLOWED_DEFINED, "user-defined alignment")
        canonical_source = {name: user[name] for name in order}
        native = dict(canonical_source)
        paml_pre = dict(canonical_source)
        coordinate_system = "user_defined_codon_alignment"
        source_path = r0 / f"00_{gene}.common_species.fasta"

    length = validate_equal_alignment(native, order, ALLOWED_DEFINED, "native canonical alignment")
    validate_equal_alignment(paml_pre, order, ALLOWED_DEFINED, "PAML pre-mask alignment")

    # Hard pre-PAML gate on the complete ORFs (see validate_complete_orf_alignment).
    # Deliberately placed before the PAML-safe view, the STOP registry and every
    # output file: nothing downstream should exist if a complete reading frame
    # was corrupted by alignment.
    orf_rows = read_tsv(r0 / f"00_{gene}.orf_status.tsv")
    complete_species = [r["species"] for r in orf_rows
                        if str(r.get("complete_orf", "")).strip().lower() == "true"]
    input_path = r0 / f"00_{gene}.common_species.gapless_for_macse.fasta"
    if not input_path.exists():
        input_path = r0 / f"00_{gene}.common_species.fasta"
    input_raw, _input_order = read_fasta(input_path)
    input_seqs = {name: seq.replace("-", "") for name, seq in input_raw.items()}
    validate_complete_orf_alignment(gene, canonical_source, native, order, complete_species,
                                    input_seqs, out, coordinate_system,
                                    on_violation=args.on_complete_orf_violation)

    if args.alignment_mode == "perform":
        # MACSE ran: join raw STOP calls to MACSE frame-phase diagnostics and
        # remap them into canonical coordinates through the four validation gates.
        raw_stops = load_raw_stop_diagnostics(r0, r1, gene)
        stop_registry = build_stop_registry(gene, raw_stops, canonical_source, user, coordinate_system)
    else:
        # --alignment defined: no MACSE. Read premature stops straight off the
        # authoritative codon alignment (see scan_defined_stops).
        stop_registry = scan_defined_stops(gene, native, order, coordinate_system)
    paml_safe, technical_stops = mask_paml_stops(paml_pre, order, stop_registry, gene, coordinate_system)
    all_stops = stop_registry + technical_stops

    # Diagnostic record of MACSE frame-restoration placeholders.  The row states
    # deliberately avoid claiming insertion/deletion direction.
    frame_rows = []
    if args.alignment_mode == "perform":
        for species in order:
            for column, ch in enumerate(macse[species], start=1):
                if ch == "!":
                    frame_rows.append({
                        "gene": gene, "species": species, "alignment_column": column,
                        "source_character": "!", "native_character": "-", "paml_safe_character": "N",
                        "interpretation": "MACSE_partial_codon_placeholder_direction_inferred_from_alignment_and_tree",
                    })

    aa = aa_alignment(paml_safe, order)
    write_fasta(native, out / f"02_{gene}.primary_codon_alignment_native.fasta", order)
    write_fasta(paml_safe, out / f"02_{gene}.primary_codon_alignment.fasta", order)
    write_fasta(paml_safe, out / f"02_{gene}.nt_alignment.fasta", order)
    write_fasta(aa, out / f"02_{gene}.primary_AA_alignment.fasta", order)
    write_phylip(paml_safe, out / f"02_{gene}.codon_for_paml.phy", order)
    Phylo.write(tree, str(out / f"02_{gene}.tree_for_asr.nwk"), "newick")

    write_tsv(all_stops, out / f"02_{gene}.masked_inframe_premature_stops_after_macse_correction.tsv", [
        "gene", "species", "codon_position", "nt_start", "nt_end",
        "primary_alignment_start", "primary_alignment_end", "stop_codon", "corrected_homologous_codon",
        "coordinate_system", "terminal_codon", "pseudogenizing_event_candidate", "independent_stop_candidate",
        "mapped_columns_contiguous", "codon_frame_aligned", "codon_confirmed_stop", "frame_shifted_at_stop", "upstream_macse_marker_count",
        "upstream_macse_frame_correction_mod3", "stop_phase_interpretation", "reason", "stop_event_key",
    ])
    write_tsv(frame_rows, out / f"02_{gene}.macse_frameshift_placeholders.tsv", [
        "gene", "species", "alignment_column", "source_character", "native_character",
        "paml_safe_character", "interpretation",
    ])
    write_tsv([{
        "gene": gene,
        "alignment_mode": args.alignment_mode,
        "aligner": "macse" if args.alignment_mode == "perform" else "NA",
        "canonical_alignment_source": str(source_path),
        "coordinate_system": coordinate_system,
        "alignment_length": length,
        "n_species": len(order),
        "defined_alignment_columns_modified": False,
        "second_alignment_performed": False,
        "premature_stop_detection": (
            "macse_raw_stop_remap_with_four_gates" if args.alignment_mode == "perform"
            else "direct_codon_frame_scan_of_defined_alignment_no_macse"
        ),
        "paml_stop_masking": "exact STOPs -> NNN; premature STOP alignment coordinates retained",
    }], out / f"02_{gene}.alignment_provenance.tsv")
    write_tsv([
        {"alignment_site": i, "codon_site": (i + 2) // 3, "coordinate_system": coordinate_system}
        for i in range(1, length + 1)
    ], out / f"02_{gene}.primary_site_map.tsv")

    n_candidate_stops = sum(1 for r in stop_registry if trueish(r.get("pseudogenizing_event_candidate")))
    print(
        f"Finished step02 for {gene}; mode={args.alignment_mode}; canonical={coordinate_system}; "
        f"species={len(order)}; columns={length}; MACSE_placeholders={len(frame_rows)}; "
        f"premature_stops_registered={len(stop_registry)}; independent_stop_candidates={n_candidate_stops}"
    )


if __name__ == "__main__":
    main()
