#!/usr/bin/env python3
"""Integrate PAML substitution ASR with Pensieve event states and reconstruct ORF history.

IndelMaP is deliberately absent from this authoritative integration.  The PAML
marginal sequence is the nucleotide scaffold for each non-root internode; the
Pensieve event engine decides structural gap/STOP states.  The dated biological
root has no distinct PAML marginal sequence under clock=0 and is therefore never
fabricated.
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

STOPS = {"TAA", "TAG", "TGA"}
ACGT = set("ACGT")


def read_fasta(path):
    records, order, name, buf = {}, [], None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(buf).upper().replace("U", "T")
                name = line[1:].split()[0]
                order.append(name); buf = []
            else:
                buf.append(line)
    if name is not None:
        records[name] = "".join(buf).upper().replace("U", "T")
    return records, order


def write_fasta(records, path, order, width=80):
    with open(path, "w") as handle:
        for name in order:
            if name not in records:
                continue
            handle.write(f">{name}\n")
            seq = records[name]
            for i in range(0, len(seq), width):
                handle.write(seq[i:i + width] + "\n")


def read_tsv(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(rows, path, header=None):
    if header is None:
        header = []
        for row in rows:
            for key in row:
                if key not in header:
                    header.append(key)
    with open(path, "w", newline="") as handle:
        w = csv.DictWriter(handle, delimiter="\t", fieldnames=header, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def safe_int(value):
    try:
        return int(float(str(value)))
    except Exception:
        return None


def trueish(value):
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


class Node:
    __slots__ = ("name", "branch_length", "children", "parent", "is_tip")
    def __init__(self):
        self.name = ""; self.branch_length = None; self.children = []; self.parent = None; self.is_tip = False


def parse_newick(path):
    text = Path(path).read_text()
    out, depth = [], 0
    for ch in text:
        if ch == "[": depth += 1
        elif ch == "]": depth = max(0, depth - 1)
        elif depth == 0: out.append(ch)
    text = "".join(out).strip().split(";")[0]
    tokens = re.findall(r"'[^']*'|\"[^\"]*\"|[(),:]|[^(),:;\s]+", text)
    pos = 0
    def unquote(x):
        if len(x) >= 2 and x[0] == x[-1] and x[0] in "'\"": return x[1:-1]
        return x
    def parse():
        nonlocal pos
        n = Node()
        if pos < len(tokens) and tokens[pos] == "(":
            pos += 1
            while True:
                c = parse(); c.parent = n; n.children.append(c)
                if pos < len(tokens) and tokens[pos] == ",": pos += 1; continue
                break
            if pos >= len(tokens) or tokens[pos] != ")": raise SystemExit(f"Malformed Newick {path}")
            pos += 1
        if pos < len(tokens) and tokens[pos] not in {",", ")", ":"}:
            n.name = unquote(tokens[pos]); pos += 1
        if pos < len(tokens) and tokens[pos] == ":":
            pos += 1; n.branch_length = float(tokens[pos]); pos += 1
        n.is_tip = not n.children
        return n
    root = parse()
    if pos != len(tokens): raise SystemExit(f"Trailing Newick content in {path}")
    return root


def iter_nodes(root, order="preorder"):
    if order == "preorder":
        stack = [root]
        while stack:
            n = stack.pop(); yield n; stack.extend(reversed(n.children))
    else:
        stack, out = [(root, False)], []
        while stack:
            n, done = stack.pop()
            if done: out.append(n)
            else:
                stack.append((n, True)); stack.extend((c, False) for c in reversed(n.children))
        yield from out


def descendant_key(node):
    return "|".join(sorted(n.name for n in iter_nodes(node) if n.is_tip))


def node_depth(node):
    depth = 0
    cur = node
    while cur.parent is not None:
        depth += 1
        cur = cur.parent
    return depth


def descendant_tips_phylo(node):
    return [n.name for n in iter_nodes(node) if n.is_tip]


def orf_check(aligned_sequence, structural_uncertain=False):
    seq = aligned_sequence.replace("-", "")
    codons = [seq[i:i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
    starts = len(seq) >= 3 and seq[:3] == "ATG"
    mod3 = len(seq) % 3 == 0
    terminal = codons[-1] if codons and mod3 else "NA"
    has_terminal_stop = terminal in STOPS
    scan = codons[:-1] if has_terminal_stop else codons
    internal_stops = [(i + 1, c) for i, c in enumerate(scan) if c in STOPS]
    n_unknown_bases = sum(ch not in ACGT for ch in seq)
    definite_disruption = (not starts) or (not mod3) or bool(internal_stops)
    sequence_uncertain = n_unknown_bases > 0
    if definite_disruption:
        status = "disrupted"
    elif structural_uncertain or sequence_uncertain:
        status = "uncertain"
    else:
        status = "intact"
    return {
        "cds_length": len(seq), "length_mod_3": len(seq) % 3,
        "starts_with_atg": starts, "has_terminal_stop": has_terminal_stop,
        "n_internal_stops": len(internal_stops),
        "internal_stop_codons": ",".join(f"{i}:{c}" for i, c in internal_stops) or "NA",
        "n_unknown_bases": n_unknown_bases,
        "structural_state_uncertain": structural_uncertain,
        "sequence_uncertain": sequence_uncertain,
        "coding_status": status,
        "coding_intact": True if status == "intact" else (False if status == "disrupted" else "NA"),
    }


def own_terminal_codon(aligned_sequence):
    """A tip's own true terminal codon, from ITS OWN sequence, ignoring the
    multiple-alignment's column coordinates entirely.

    A tip's real, biological sequence has no gap in it at all; the gap only
    exists in the shared multiple-alignment coordinate system, inserted
    wherever OTHER lineages have length this tip does not. When another
    lineage has unrelated extra sequence anywhere near a tip's own true
    C-terminus, that padding can land in the middle of -- or immediately
    after -- this tip's own last codon, scattering its 3 real bases across
    non-adjacent alignment columns even though, in the tip's own real
    sequence, they are perfectly contiguous. Reading a fixed span of raw
    alignment columns is therefore not reliable for finding "the" terminal
    codon; stripping gaps first and taking the last 3 real bases is.
    Returns a 3-character codon, or None if fewer than 3 real bases exist.
    """
    ungapped = aligned_sequence.replace("-", "")
    return ungapped[-3:].upper() if len(ungapped) >= 3 else None


def fitch_reconstruct_codon(root, tip_codons):
    """Two-pass Fitch parsimony treating a whole codon as one discrete state.

    PAML's codon-substitution model has a 61-state (sense-codon-only) state
    space by construction, so it can never reconstruct a stop codon at any
    node -- including the ordinary, universal terminal stop present in
    virtually every intact tip, which is not a lesion at all. That column is
    masked to N for every tip before codeml ever sees it, so PAML's own
    reconstruction there is not just wrong, it is uninformed: nothing in its
    model or its input constrains it. This reconstructs that span directly
    from the real observed tip data instead.

    Ordinary bottom-up parsimony: look at the state of the tips (or already-
    resolved children) under a node and assign that node the same state when
    they agree, all the way from the tips up to the root. The state used per
    tip is its own derived terminal codon (see own_terminal_codon()) rather
    than a fixed span of raw alignment columns, so a real, contiguous codon
    that a nearby unrelated insertion happened to scatter across non-adjacent
    columns is still read correctly.

    tip_codons: {tip_name: 3-character codon string, or None if undetermined}.
    Returns {node_name: resolved_codon_or_None} for every node in the tree.
    """
    down = {}
    for node in iter_nodes(root, order="postorder"):
        if node.is_tip:
            codon = tip_codons.get(node.name)
            down[node] = {codon} if codon else None
        else:
            child_sets = [s for s in (down[c] for c in node.children) if s is not None]
            if not child_sets:
                down[node] = None
            else:
                inter = set.intersection(*child_sets)
                down[node] = inter if inter else set.union(*child_sets)
    resolved = {}
    for node in iter_nodes(root, order="preorder"):
        candidates = down[node]
        parent_codon = resolved.get(node.parent.name) if node.parent is not None else None
        if candidates is None:
            resolved[node.name] = None
        elif parent_codon is not None and parent_codon in candidates:
            resolved[node.name] = parent_codon
        else:
            resolved[node.name] = min(candidates)
    return resolved


def main():
    parser = argparse.ArgumentParser(description="Build lesion-aware ancestral sequences and ORF history.")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--paml-tree", required=True)
    parser.add_argument("--pensieve-tree", required=True)
    parser.add_argument("--ancestral-fasta", required=True, help="PAML marginal internal-node FASTA")
    parser.add_argument("--tip-alignment", required=True, help="Canonical native alignment for observed tips")
    parser.add_argument("--character-states", required=True)
    parser.add_argument("--characters", required=True)
    parser.add_argument("--events", required=True)
    parser.add_argument("--frame-arithmetic", default=None)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    gene = args.gene
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    pensieve = parse_newick(args.pensieve_tree)
    paml_tree = parse_newick(args.paml_tree)
    paml_asr, _ = read_fasta(args.ancestral_fasta)
    tips_native, _ = read_fasta(args.tip_alignment)

    pens_by_key = {descendant_key(n): n.name for n in iter_nodes(pensieve) if not n.is_tip}
    paml_by_key = {descendant_key(n): n.name for n in iter_nodes(paml_tree) if not n.is_tip and n.parent is not None}
    rename = {}
    for key, plabel in paml_by_key.items():
        if key not in pens_by_key:
            raise SystemExit(f"PAML-labelled reporting node {plabel} has no Pensieve descendant-set match: {key}")
        rename[plabel] = pens_by_key[key]
    write_tsv([
        {"gene": gene, "paml_node_label": p, "pensieve_node_label": q,
         "rooted_descendant_tip_key": key}
        for key, p in sorted(paml_by_key.items()) for q in [rename[p]]
    ], out / f"03_{gene}.paml_to_pensieve_node_map.tsv",
       ["gene", "paml_node_label", "pensieve_node_label", "rooted_descendant_tip_key"])

    scaffold = {rename[p]: seq for p, seq in paml_asr.items() if p in rename}
    node_order = [n.name for n in iter_nodes(pensieve)]
    root_name = pensieve.name
    if root_name in scaffold:
        scaffold.pop(root_name, None)  # never accept an accidental fabricated root sequence

    characters = {r["character_id"]: r for r in read_tsv(args.characters)}
    states = {(r["character_id"], r["node_label"]): r.get("state", "unknown") for r in read_tsv(args.character_states)}
    events = read_tsv(args.events)

    # Reconstruct the terminal (biologically ordinary, non-disrupting) stop
    # codon at internal nodes directly from each tip's own true sequence,
    # once, up front -- see own_terminal_codon()/fitch_reconstruct_codon()
    # for why neither PAML's own scaffold nor a fixed span of raw alignment
    # columns can be trusted for this. The result is written into the fixed
    # final 3 columns of the "integrated" alignment below purely to keep
    # every record the same width; the VALUE itself comes from each
    # descendant's own real sequence, not from those columns.
    tip_align_lengths = {len(seq) for seq in tips_native.values()}
    terminal_fitch_columns = []
    if len(tip_align_lengths) == 1:
        aln_len = next(iter(tip_align_lengths))
        term_start0 = aln_len - 3
        if term_start0 >= 0:
            tip_codons = {t: own_terminal_codon(seq) for t, seq in tips_native.items()}
            resolved_codons = fitch_reconstruct_codon(pensieve, tip_codons)
            for offset in range(3):
                col = term_start0 + offset
                resolved_bases = {label: codon[offset] for label, codon in resolved_codons.items() if codon}
                terminal_fitch_columns.append((col, resolved_bases))

    integrated = {}
    native_cds = {}
    status_rows = []
    status_by_node = {}
    structural_uncertain_by_node = {}

    for node in iter_nodes(pensieve):
        label = node.name
        if node.is_tip:
            seq = tips_native.get(label)
            if seq is None:
                raise SystemExit(f"Tip {label} missing from native alignment")
            integrated[label] = seq
            uncertain = False
        else:
            seq = scaffold.get(label)
            if seq is None:
                status_by_node[label] = "unavailable"
                structural_uncertain_by_node[label] = True
                status_rows.append({
                    "gene": gene, "node_label": label, "node_type": "root" if node.parent is None else "internal",
                    "sequence_source": "no_distinct_paml_clock0_root_sequence" if node.parent is None else "missing_paml_sequence",
                    "cds_length": "NA", "length_mod_3": "NA", "starts_with_atg": "NA",
                    "has_terminal_stop": "NA", "n_internal_stops": "NA", "internal_stop_codons": "NA",
                    "n_unknown_bases": "NA", "structural_state_uncertain": True,
                    "sequence_uncertain": True, "coding_status": "unavailable", "coding_intact": "NA",
                })
                continue
            chars = list(seq)
            uncertain = False
            for cid, char in characters.items():
                state = states.get((cid, label), "unknown")
                start = safe_int(char.get("alignment_start")); end = safe_int(char.get("alignment_end"))
                if start is None or end is None:
                    continue
                if state in {"ambiguous", "unknown", "NA", ""}:
                    # An ambiguous/unresolved reconstruction only makes this
                    # node's ORF-completeness call itself uncertain when the
                    # character could actually change that call: a STOP
                    # allele (is a stop present or not), or an indel whose
                    # length shifts the reading frame. A gap-vs-residue tie
                    # on an in-frame indel changes neither the frame nor the
                    # stop-codon content of the resulting sequence, so it
                    # cannot flip intact/disrupted; treating every such tie
                    # as grounds to mark the whole node "uncertain" is what
                    # made genes with many small in-frame indels (a normal,
                    # non-disrupting occurrence) read as almost entirely
                    # amber/uncertain even where every actual ORF check below
                    # is clean.
                    klass = char.get("character_class")
                    length_mod3 = safe_int(char.get("length_mod_3"))
                    orf_relevant = klass == "stop_mask" or (klass == "indel" and length_mod3 not in {None, 0})
                    if orf_relevant:
                        uncertain = True
                    continue
                if char.get("character_class") == "indel":
                    if state == "gap":
                        for i in range(start - 1, min(end, len(chars))): chars[i] = "-"
                elif char.get("character_class") == "stop_mask" and state == "stop_present":
                    codon = str(char.get("stop_codon", "")).upper()
                    if codon not in STOPS:
                        uncertain = True
                    else:
                        fill = (codon * (((end - start + 1) // 3) + 1))[:end - start + 1]
                        for offset, ch in enumerate(fill):
                            if start - 1 + offset < len(chars): chars[start - 1 + offset] = ch
            for col, resolved in terminal_fitch_columns:
                if col < len(chars) and label in resolved:
                    chars[col] = resolved[label]
            integrated[label] = "".join(chars)

        result = orf_check(integrated[label], structural_uncertain=uncertain)
        native_cds[label] = integrated[label].replace("-", "")
        status_by_node[label] = result["coding_status"]
        structural_uncertain_by_node[label] = uncertain
        status_rows.append({
            "gene": gene, "node_label": label,
            "node_type": "tip" if node.is_tip else ("root" if node.parent is None else "internal"),
            "sequence_source": "observed_native_tip" if node.is_tip else "paml_scaffold_plus_pensieve_event_states",
            **result,
        })

    # All inspectable ancestral/tip sequence outputs follow the rooted tree in
    # deterministic depth-first preorder: an internode is immediately followed
    # by the clade descending from it, and sister clades follow the left-to-right
    # order of the user-supplied tree. Missing root sequence under clock=0 is
    # simply skipped by write_fasta without disturbing the remaining order.
    write_fasta(scaffold, out / f"03_{gene}.paml_scaffold_relabelled.fa", node_order)
    write_fasta(integrated, out / f"03_{gene}.ancestral_integrated_alignment.fa", node_order)
    write_fasta(integrated, out / f"03_{gene}.phylogenetic_msa.fa", node_order)
    write_fasta(integrated, out / f"03_{gene}.ancestral_native.fa", node_order)  # compatibility alias
    write_fasta(native_cds, out / f"03_{gene}.ancestral_native_cds.fa", node_order)

    order_rows = []
    sequence_rank = 0
    for phylo_rank, node in enumerate(iter_nodes(pensieve), start=1):
        available = node.name in integrated
        if available:
            sequence_rank += 1
        desc = descendant_tips_phylo(node)
        order_rows.append({
            "gene": gene,
            "phylogenetic_rank": phylo_rank,
            "sequence_fasta_rank": sequence_rank if available else "NA",
            "node_label": node.name,
            "node_type": "tip" if node.is_tip else ("root" if node.parent is None else "internal"),
            "depth_from_root": node_depth(node),
            "parent_label": node.parent.name if node.parent is not None else "NA",
            "sequence_available": available,
            "descendant_tip_count": len(desc),
            "descendant_tips_in_phylogenetic_order": ",".join(desc),
        })
    write_tsv(order_rows, out / f"03_{gene}.phylogenetic_sequence_order.tsv", [
        "gene", "phylogenetic_rank", "sequence_fasta_rank", "node_label", "node_type",
        "depth_from_root", "parent_label", "sequence_available", "descendant_tip_count",
        "descendant_tips_in_phylogenetic_order",
    ])
    write_tsv(status_rows, out / f"04_{gene}.ancestral_orf_walk.tsv")

    events_by_branch = defaultdict(list)
    for ev in events:
        events_by_branch[ev.get("branch", "")].append(ev)

    # A disabling state may already be present at the sampled root.  In that
    # case its origin predates the tree: descendants inherit pseudogenic history
    # but no within-tree branch is falsely labelled as the first loss.  Root
    # ambiguity is carried as unknown rather than guessed away.
    root_disabling = False
    root_history_uncertain = False
    for cid, char in characters.items():
        root_state = states.get((cid, pensieve.name), "unknown")
        if root_state in {"ambiguous", "unknown", "NA", ""}:
            if char.get("character_class") == "stop_mask":
                root_history_uncertain = True
            elif char.get("character_class") == "indel":
                length = safe_int(char.get("length_mod_3"))
                if length not in {None, 0}:
                    root_history_uncertain = True
            continue
        if char.get("character_class") == "indel":
            length = safe_int(char.get("length_mod_3"))
            if root_state == "gap" and length not in {None, 0}:
                # A root GAP is an ancestral occupancy state, not evidence that
                # a deletion happened before the sampled root: the alternative
                # history may be an insertion on the residue lineage.  Without
                # a distinct clock=0 PAML root sequence we therefore mark the
                # pre-root pseudogenic history as unresolved, never as known.
                root_history_uncertain = True
        elif char.get("character_class") == "stop_mask" and root_state == "stop_present":
            # An independently mapped premature STOP allele is itself a
            # disabling nucleotide state, so confident root presence supports
            # a loss predating the sampled tree.
            root_disabling = True

    known_history = {pensieve.name: True if root_disabling else (None if root_history_uncertain else False)}
    branch_rows = []
    for node in iter_nodes(pensieve):
        if node.parent is None:
            continue
        parent = node.parent.name; child = node.name; branch = f"{parent}->{child}"
        parent_status = status_by_node.get(parent, "unavailable")
        child_status = status_by_node.get(child, "unavailable")
        acquired = events_by_branch.get(branch, [])
        confident_disabling = []
        ambiguous_disabling = []
        for ev in acquired:
            is_stop_gain = ev.get("character_class") == "stop_mask" and ev.get("biological_interpretation") == "premature_stop_gained"
            is_frameshift = ev.get("character_class") == "indel" and safe_int(ev.get("length_mod_3")) not in {None, 0}
            if not (is_stop_gain or is_frameshift):
                continue
            if trueish(ev.get("direction_confident")):
                confident_disabling.append(ev.get("event_id", "NA"))
            else:
                ambiguous_disabling.append(ev.get("event_id", "NA"))

        parent_hist = known_history.get(parent, None)
        root_adjacent = node.parent is pensieve
        uncertainty_reason = "NA"
        transition_evidence = "NA"

        if parent_hist is True:
            # Once a disabling history is known upstream it remains sticky. A
            # compensating lesion may restore frame/apparent ORF but does not
            # resurrect the ancestral protein or create another first-loss call.
            child_hist = True
            if child_status == "intact":
                transition = "apparent_orf_restoration"
                transition_evidence = "inherited_pseudogenic_history_plus_current_orf_intact"
            else:
                transition = "already_pseudogenic"
                transition_evidence = "inherited_pseudogenic_history"

        elif parent_hist is None:
            # Entering pseudogenic history is unresolved. A definite lesion can
            # prove the child has a disabling history, but cannot prove this is
            # the FIRST loss on the sampled tree.
            if confident_disabling:
                child_hist = True
                transition = "confirmed_disabling_event_first_loss_unresolved"
                transition_evidence = "confident_disabling_event"
                uncertainty_reason = "entering_pseudogenic_history_unresolved"
            elif child_status == "disrupted":
                child_hist = True
                transition_evidence = "child_sequence_disrupted_without_confident_catalogued_branch_event"
                if root_adjacent:
                    transition = "root_adjacent_disruption_first_loss_unresolved"
                    uncertainty_reason = "no_distinct_clock0_biological_root_sequence"
                else:
                    transition = "sequence_disruption_first_loss_unresolved"
                    uncertainty_reason = "entering_pseudogenic_history_unresolved"
            elif ambiguous_disabling:
                child_hist = None
                transition = "ambiguous_disabling_event"
                transition_evidence = "ambiguous_disabling_event"
                uncertainty_reason = "event_origin_or_direction_ambiguous"
            elif child_status == "intact":
                # This node's own reconstructed sequence is confidently intact
                # and no disabling event (confident or ambiguous) was
                # catalogued on this exact branch. Whatever remained
                # unresolved about pre-tree/root-adjacent history, this node
                # is direct, confident evidence of an intact ORF at this
                # point, so the inherited uncertainty is resolved to "clean"
                # from here rather than propagated indefinitely down every
                # descendant branch regardless of their own, individually
                # confident, intact status.
                child_hist = False
                transition = "intact"
                transition_evidence = "no_inferred_disabling_history"
            else:
                child_hist = None
                transition = "history_uncertain_no_confirmed_event"
                transition_evidence = "no_confident_disabling_event"
                uncertainty_reason = "entering_pseudogenic_history_unresolved"

        elif confident_disabling:
            child_hist = True
            transition = "pseudogenization"
            transition_evidence = "confident_disabling_event"

        elif parent_status == "intact" and child_status == "disrupted":
            # Captures sequence-level disabling changes not represented by the
            # indel/STOP catalogue (for example a reconstructed start-codon loss).
            child_hist = True
            transition = "pseudogenization"
            transition_evidence = "parent_intact_child_disrupted_sequence_state"

        elif ambiguous_disabling:
            child_hist = None
            transition = "ambiguous_disabling_event"
            transition_evidence = "ambiguous_disabling_event"
            uncertainty_reason = "event_origin_or_direction_ambiguous"

        elif child_status == "disrupted":
            # Most importantly, a branch immediately below the clock=0 biological
            # root cannot be claimed as the first substitution-only loss because
            # no distinct root sequence exists.  The child disruption is real,
            # however, so descendants inherit known pseudogenic history.
            child_hist = True
            transition_evidence = "child_sequence_disrupted_without_confident_catalogued_branch_event"
            if root_adjacent:
                transition = "root_adjacent_disruption_first_loss_unresolved"
                uncertainty_reason = "no_distinct_clock0_biological_root_sequence"
            else:
                transition = "sequence_disruption_first_loss_unresolved"
                uncertainty_reason = "parent_sequence_not_confidently_intact"

        elif child_status == "uncertain":
            child_hist = None
            transition = "sequence_state_uncertain"
            transition_evidence = "child_sequence_uncertain"
            uncertainty_reason = "ancestral_sequence_or_structural_state_uncertain"

        elif child_status == "unavailable":
            child_hist = None
            transition = "sequence_state_unavailable"
            transition_evidence = "child_sequence_unavailable"
            uncertainty_reason = "no_reconstructable_child_sequence"

        else:
            child_hist = False
            transition = "intact"
            transition_evidence = "no_inferred_disabling_history"

        known_history[child] = child_hist

        branch_rows.append({
            "gene": gene, "branch": branch, "parent_node": parent, "child_node": child,
            "child_is_tip": node.is_tip, "parent_coding_status": parent_status,
            "child_coding_status": child_status, "orf_transition": transition,
            "transition_evidence": transition_evidence,
            "uncertainty_reason": uncertainty_reason,
            "known_pseudogenic_history": "unknown" if child_hist is None else child_hist,
            "confident_disabling_events_on_branch": ",".join(confident_disabling) or "NA",
            "ambiguous_disabling_events_on_branch": ",".join(ambiguous_disabling) or "NA",
            "all_events_on_branch": ",".join(ev.get("event_id", "NA") for ev in acquired) or "NA",
            "n_events_on_branch": len(acquired),
        })
    write_tsv(branch_rows, out / f"04_{gene}.orf_transitions_by_branch.tsv")

    # Frame cross-check compares current frame offset with reconstructed CDS
    # length modulo three only. It intentionally does not compare frame state to
    # overall ORF status because nonsense mutations do not change frame.
    mismatches = []
    if args.frame_arithmetic:
        frame = {r.get("node_label"): r for r in read_tsv(args.frame_arithmetic)}
        by_status = {r["node_label"]: r for r in status_rows}
        for label, fr in frame.items():
            sr = by_status.get(label)
            if not sr or sr.get("length_mod_3") in {"NA", None, ""}:
                continue
            if trueish(fr.get("structural_state_ambiguous")):
                continue
            expected_shift = safe_int(sr.get("length_mod_3")) != 0
            frame_shift = trueish(fr.get("frame_currently_shifted", fr.get("frame_disrupted_provisional")))
            if expected_shift != frame_shift:
                mismatches.append({
                    "gene": gene, "node_label": label, "frame_arithmetic_shifted": frame_shift,
                    "integrated_cds_length_mod_3": sr.get("length_mod_3"),
                    "note": "signed frame arithmetic and integrated native CDS length disagree",
                })
    write_tsv(mismatches, out / f"04_{gene}.orf_frame_crosscheck_mismatches.tsv",
              ["gene", "node_label", "frame_arithmetic_shifted", "integrated_cds_length_mod_3", "note"])

    print(f"Finished ORF/history integration for {gene}; nodes_with_sequence={len(integrated)}; "
          f"pseudogenization_branches={sum(r['orf_transition']=='pseudogenization' for r in branch_rows)}; "
          f"apparent_restorations={sum(r['orf_transition']=='apparent_orf_restoration' for r in branch_rows)}")


if __name__ == "__main__":
    main()
