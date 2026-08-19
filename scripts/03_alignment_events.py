#!/usr/bin/env python3
"""Breakpoint-coded indel and mask events, placed on the dated tree by parsimony.

Design notes (v3.30)
--------------------
1. INDEL CHARACTERS ARE CODED BY BREAKPOINTS, NOT BY COLUMNS.
   Earlier versions cut a new block whenever the exact carrier set changed
   between adjacent columns, so any unrelated lineage with an overlapping indel
   split someone else's event.  On GUCA1B that fragmented one 42 bp Miniopterus
   deletion into three blocks, and it can turn a single in-frame deletion into
   several apparent frameshifts.  Here each maximal gap run per tip is extracted
   as an interval, intervals are clustered by shared breakpoints, and each
   cluster is ONE binary character (Simmons & Ochoterena simple indel coding).

2. MACSE '!' IS A PARTIAL-CODON PLACEHOLDER, NOT AN EVENT DIRECTION.
   In the native structural view '!' is rendered as '-' so the artificial frame-
   restoration placeholder is removed. Whether the homologous change is an
   insertion or deletion is inferred from the complete aligned occupancy pattern
   and the phylogeny; Pensieve never equates one '!' with one deleted nucleotide.

3. PLACEMENT IS BY PARSIMONY, NOT BY MRCA OF CARRIERS.
   Equal-cost binary Sankoff on the dated tree.  Origins are branch transitions.
   Root polarity is reconstructed, so a state carried by almost every tip is
   reported as ancestral with descendant losses rather than as a root event.

4. SUPPORT IS REPORTED.
   For every node the cost of forcing PRESENT vs ABSENT is recorded; the
   difference (delta-parsimony) is how many extra steps the alternative history
   costs.  delta == 0 is an exact tie and is flagged, never silently resolved.

5. FRAME ARITHMETIC IS SIGNED AND HISTORY-AWARE.
   Confident residue->gap transitions contribute negative length (deletion) and
   gap->residue transitions positive length (insertion). Current frame offset and
   prior frameshifting history are reported separately. STOP mutations are not
   called frame disruptions.

Standard library only: no Biopython, no network.  This keeps the event engine
testable in isolation from codeml.
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

ABSENT, PRESENT = 0, 1
INF = float("inf")
ACGT = set("ACGT")


# ---------------------------------------------------------------- io

def read_fasta(path):
    records, order, name, buf = {}, [], None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(buf)
                header = line[1:].strip()
                if not header:
                    raise SystemExit(f"{path}: FASTA record with an empty identifier")
                name = header.split()[0]
                if name in records:
                    raise SystemExit(f"{path}: duplicate FASTA identifier {name!r}")
                order.append(name)
                buf = []
            elif line.strip():
                buf.append(line.strip())
    if name is not None:
        records[name] = "".join(buf)
    if not records:
        raise SystemExit(f"{path}: no sequences found")
    records = {k: v.upper().replace("U", "T") for k, v in records.items()}
    lengths = {len(v) for v in records.values()}
    if len(lengths) != 1:
        raise SystemExit(f"{path}: sequences are not all the same length {sorted(lengths)}")
    return records, order, lengths.pop()


def read_tsv(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(rows, path, header=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if header is None:
        header = []
        for row in rows:
            for key in row:
                if key not in header:
                    header.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(header),
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def safe_int(value):
    try:
        return int(float(str(value)))
    except Exception:
        return None


def trueish(value):
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


# ---------------------------------------------------------------- tree

class Node:
    __slots__ = ("name", "branch_length", "children", "parent", "is_tip", "idx")

    def __init__(self):
        self.name = ""
        self.branch_length = None
        self.children = []
        self.parent = None
        self.is_tip = False
        self.idx = -1


def parse_newick(path):
    text = Path(path).read_text()
    out, depth = [], 0
    for ch in text:                      # strip [comments]
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth = max(0, depth - 1)
        elif depth == 0:
            out.append(ch)
    text = "".join(out).strip()
    if ";" in text:
        text = text[:text.index(";")]
    tokens = re.findall(r"'[^']*'|\"[^\"]*\"|[(),:]|[^(),:;\s]+", text.strip())
    pos = 0

    def unquote(label):
        label = label.strip()
        if len(label) >= 2 and label[0] == label[-1] and label[0] in "'\"":
            label = label[1:-1]
        return label.strip()

    def parse():
        nonlocal pos
        node = Node()
        if pos < len(tokens) and tokens[pos] == "(":
            pos += 1
            while True:
                child = parse()
                child.parent = node
                node.children.append(child)
                if pos < len(tokens) and tokens[pos] == ",":
                    pos += 1
                    continue
                break
            if pos >= len(tokens) or tokens[pos] != ")":
                raise SystemExit(f"{path}: malformed Newick, expected ')'")
            pos += 1
        if pos < len(tokens) and tokens[pos] not in {",", ")", ":"}:
            node.name = unquote(tokens[pos])
            pos += 1
        if pos < len(tokens) and tokens[pos] == ":":
            pos += 1
            try:
                node.branch_length = float(tokens[pos])
            except (IndexError, ValueError):
                raise SystemExit(f"{path}: malformed branch length")
            pos += 1
        node.is_tip = not node.children
        return node

    root = parse()
    if pos != len(tokens):
        raise SystemExit(f"{path}: trailing content after the root clade")
    return root


def iter_nodes(root, order="preorder"):
    if order == "preorder":
        stack = [root]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))
    else:
        out, stack = [], [(root, False)]
        while stack:
            node, done = stack.pop()
            if done:
                out.append(node)
            else:
                stack.append((node, True))
                stack.extend((c, False) for c in reversed(node.children))
        yield from out


def collapse_unifurcations(root):
    removed = []
    changed = True
    while changed:
        changed = False
        for node in list(iter_nodes(root, "preorder")):
            if node.is_tip or len(node.children) != 1:
                continue
            child = node.children[0]
            child.branch_length = (node.branch_length or 0.0) + (child.branch_length or 0.0)
            child.parent = node.parent
            if node.parent is None:
                root = child
            else:
                node.parent.children[node.parent.children.index(node)] = child
            removed.append(node.name or "<unlabelled>")
            changed = True
            break
    return root, removed


def apply_pensieve_labels(root, prefix="Node"):
    """Label EVERY internal node with Pensieve's own postorder numbering.

    Pensieve owns the internode namespace.  codeml's numbers are mapped onto
    these labels by descendant-tip set in step 03b, so labels stay stable and
    defined even when codeml is not run.
    """
    original = {}
    counter = 0
    for node in iter_nodes(root, "postorder"):
        if node.is_tip:
            if not node.name:
                raise SystemExit("Tree contains an unnamed tip")
            continue
        counter += 1
        original[f"{prefix}{counter}"] = node.name or "NA"
        node.name = f"{prefix}{counter}"
    for i, node in enumerate(iter_nodes(root, "postorder")):
        node.idx = i
    return root, original


def write_newick(root, path):
    def render(node):
        text = f"({','.join(render(c) for c in node.children)}){node.name}" if node.children else node.name
        if node.branch_length is not None:
            text += f":{node.branch_length:.10g}"
        return text
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(render(root) + ";\n")


def node_ages(root, dated):
    """Age before present for every node, from root-to-node depth."""
    depth = {root: 0.0}
    for node in iter_nodes(root, "preorder"):
        for child in node.children:
            depth[child] = depth[node] + (child.branch_length or 0.0)
    tip_depth = [depth[n] for n in iter_nodes(root, "preorder") if n.is_tip]
    tallest = max(tip_depth) if tip_depth else 0.0
    if not dated:
        return {n.name: "NA" for n in iter_nodes(root, "preorder")}
    return {n.name: round(tallest - depth[n], 6) for n in iter_nodes(root, "preorder")}


# ---------------------------------------------------------------- parsimony

def clamp_boundary_evidence(pair, cap):
    """Bound how much a subtree's own internal cost pattern can influence
    its PARENT's Sankoff decision, regardless of how many tips the subtree
    contains.

    Real problem (see docstring on pseudogenic_components() below): a
    monophyletic pseudogenized clade with many descendant tips does not
    inflate a SINGLE character's own cost trade-off by tip count alone --
    ordinary Sankoff already handles a perfectly homogeneous clade (one
    shared gain, zero extra cost regardless of size) correctly. The real
    distortion is statistical: with more descendant tips, there are simply
    more chances for SOME subset of them to happen to share a given
    indel/STOP character densely enough that "call the clade's ancestor a
    carrier" becomes the cheaper explanation for THAT one character -- and
    with many structural characters scored independently, enough of them
    tip that way to make an ancestor look multiply-disrupted in aggregate,
    even when each individual call was locally "correct" under ordinary
    parsimony. Capping the cost DIFFERENCE a component can hand to its own
    parent (rather than the raw, unbounded pair) means a 20-tip pseudogenic
    radiation can carry at most `cap` worth of extra weight toward "this
    character was already present at the component's entry" -- the same
    ceiling a 2-3 tip version of the identical clade would already produce
    on its own, satisfying the sampling-density-invariance requirement
    without touching the component's own, fully-resolved internal history.
    """
    d0, d1 = pair
    diff = d1 - d0
    if diff > cap:
        diff = cap
    elif diff < -cap:
        diff = -cap
    return (0.0, diff) if diff >= 0 else (-diff, 0.0)


def sankoff(nodes_post, tip_states, gain_cost, loss_cost, tie_break,
            component_entry_idx=None, boundary_cap=None, pinned=None):
    """Binary Sankoff with explicit tie reporting and a representative history.

    pinned (optional; default None): a {node_idx: state} map of internal nodes
    whose ancestral state is FIXED as a hard biological constraint before
    parsimony resolves the rest of the tree. Used by the functional-consensus
    rule (see functional_consensus_pin()): when two or more phylogenetically
    independent complete-ORF (functional) lineages share exactly the same indel,
    that indel must have been present in their common functional ancestor rather
    than gained convergently in each lineage, so the character's state at that
    ancestor is pinned PRESENT. A pinned node is never reported as ambiguous.

    A representative assignment is needed for plotting/output, but equal-cost
    alternatives remain explicitly ambiguous. No tie-break setting is allowed to
    convert an exact tie into biological certainty.

    component_entry_idx/boundary_cap (both optional; default None/None
    preserves the exact original, uncapped behavior): node indices that are
    the entry point of a Stage-B pseudogenic component (see
    pseudogenic_components()). Each such node's OWN down-cost pair is still
    computed normally, in full, from its own real children -- nothing about
    the component's internal resolution changes. Only the copy of that pair
    HANDED TO THE NODE'S PARENT is bounded (clamp_boundary_evidence), so
    everything ABOVE the component boundary sees at most `boundary_cap`
    worth of evidence from the whole component, however many tips it
    actually contains.

    Capping is intentionally NOT applied in the traceback/up-pass: the
    entry node's OWN reported state and confidence use its full, unbounded
    down-cost pair, exactly like every node strictly inside the component.
    An earlier version of this correction also capped the entry node's own
    up-pass decision, which is stronger than the sampling-density fix
    requires and can overwrite real internal evidence at the entry node
    itself (e.g. GRK7/GUCY2F review). Only the MESSAGE that crosses the
    component boundary upward -- i.e. what an ancestor ABOVE the entry node
    uses when computing its own cost -- is bounded.
    """
    trans = ((0.0, gain_cost), (loss_cost, 0.0))
    down = {}
    entry_idx = component_entry_idx or set()
    pins = pinned or {}
    for node in nodes_post:
        if node.is_tip:
            obs = tip_states.get(node.name)
            down[node.idx] = [0.0, 0.0] if obs is None else (
                [INF, 0.0] if obs == PRESENT else [0.0, INF])
            continue
        costs = [0.0, 0.0]
        for parent_state in (ABSENT, PRESENT):
            costs[parent_state] = sum(
                min(trans[parent_state][s] +
                    (clamp_boundary_evidence(down[c.idx], boundary_cap)[s]
                     if c.idx in entry_idx and boundary_cap is not None else down[c.idx][s])
                    for s in (ABSENT, PRESENT))
                for c in node.children)
        # Hard constraint: a pinned internal node's own state is fixed, so the
        # subtree cost for the other state is made prohibitive (its real cost
        # given the pinned state is retained for the parent's accounting).
        if node.idx in pins:
            costs[1 - pins[node.idx]] = INF
        down[node.idx] = costs

    root = nodes_post[-1]
    rc = down[root.idx]
    best = min(rc)
    ambiguous = set()
    if root.idx in pins:
        root_state = pins[root.idx]          # fixed constraint, never ambiguous
    elif rc[ABSENT] == rc[PRESENT]:
        ambiguous.add(root.idx)
        # Representative only. The root remains 'ambiguous' in all biological
        # output fields. 'terminal' prefers the state that avoids an ancestral
        # gap; 'ancestral' prefers PRESENT; 'none' uses ABSENT deterministically.
        root_state = PRESENT if tie_break == 'ancestral' else ABSENT
    else:
        root_state = ABSENT if rc[ABSENT] < rc[PRESENT] else PRESENT

    assign = {root.idx: root_state}
    delta = {root.idx: abs(rc[PRESENT] - rc[ABSENT])}
    for node in reversed(nodes_post):
        if node.is_tip:
            continue
        pstate = assign[node.idx]
        for child in node.children:
            # The traceback/up-pass always uses each child's real, unbounded
            # down-cost pair -- including when child is itself a component
            # entry node. Bounding only applies to the down-pass sum used by
            # a node's PARENT (see sankoff()'s docstring): the entry node's
            # own C=ABSENT-vs-PRESENT decision must reflect its full
            # descendant evidence, exactly like every node strictly inside
            # the component.
            child_down = down[child.idx]
            options = [(trans[pstate][s] + child_down[s], s) for s in (ABSENT, PRESENT)]
            lo = min(v for v, _ in options)
            tied = [s for v, s in options if v == lo]
            delta[child.idx] = abs(options[0][0] - options[1][0])
            if len(tied) > 1:
                ambiguous.add(child.idx)
                if child.is_tip:
                    chosen = pstate           # missing-data tip: representative inheritance
                elif tie_break == "ancestral":
                    chosen = PRESENT
                elif tie_break == "terminal":
                    chosen = pstate
                else:
                    chosen = pstate
            else:
                chosen = tied[0]
            assign[child.idx] = chosen
    return assign, best, ambiguous, delta


# ---------------------------------------------------------- ORF-aware boundary evidence

def read_tip_orf_states(orf_status_rows, tips):
    """Definite per-tip Functional/Pseudogenic calls from 00_<gene>.orf_status
    .tsv's own complete_orf column -- the SAME per-species, reference-free
    classification already used everywhere else in Pensieve (starts with
    ATG, length a multiple of 3, no internal in-frame STOP; a MISSING
    terminal STOP is deliberately not disqualifying). A tip absent from the
    status table, or otherwise not resolvable, contributes no forced state
    (None), exactly like an unobserved tip in any other Sankoff character
    here -- it must not vote either way.
    """
    by_species = {r.get("species"): r for r in orf_status_rows}
    states = {}
    for tip in tips:
        row = by_species.get(tip)
        if row is None:
            states[tip] = None
            continue
        states[tip] = PRESENT if not trueish(row.get("complete_orf")) else ABSENT
    return states


def pseudogenic_components(nodes_post, root, tip_orf_states, orf_loss_cost, orf_restoration_cost, tie_break):
    """Stage A + Stage B: a provisional, tree-wide ORF-history map, and the
    set of branches that are the confirmed entry point of a contiguous
    pseudogenic (P) component.

    This is deliberately computed ONCE, from tip-level ORF completeness
    evidence ONLY (never from any single structural character's own
    reconstruction), specifically to avoid circularity: using one
    potentially-biased indel/STOP reconstruction to decide what counts as
    "pseudogenic," and then using that same verdict to fix the
    reconstruction it came from. F->P ("the gene broke") costs
    orf_loss_cost; the reverse, P->F ("an apparently dead gene came back"),
    costs the much larger orf_restoration_cost by default -- biologically
    rare and should not be invented for free just because a component
    happens to contain many easy, independent F->P losses instead. The root
    is never forced: sankoff() only ever reports a representative value at
    a genuine tie, and the tie itself is preserved in `ambiguous`.

    A branch parent->child is a confirmed component entry only when BOTH
    the parent's and the child's own Stage-A state are unambiguous (parent
    confidently F, child confidently P) -- an uncertain entering history is
    never upgraded into an invented first-loss boundary; it simply does not
    trigger any Stage-C compression, matching Stage B's "carry the
    uncertainty forward" requirement.

    Returns (orf_assign, orf_ambiguous, component_entry_idx, component_of)
    where component_of maps every node index inside a component (the entry
    node and everything below it, down to the next nested entry if any) to
    the entry node's own name, for audit/reporting purposes only.
    """
    orf_assign, _score, orf_ambiguous, _delta = sankoff(
        nodes_post, tip_orf_states, orf_loss_cost, orf_restoration_cost, tie_break)

    component_entry_idx = set()
    for node in iter_nodes(root, "preorder"):
        if node.parent is None:
            continue
        if node.parent.idx in orf_ambiguous or node.idx in orf_ambiguous:
            continue
        if orf_assign[node.parent.idx] == ABSENT and orf_assign[node.idx] == PRESENT:
            component_entry_idx.add(node.idx)

    component_of = {}
    for entry_node in (n for n in iter_nodes(root, "preorder") if n.idx in component_entry_idx):
        for descendant in iter_nodes(entry_node, "preorder"):
            component_of.setdefault(descendant.idx, entry_node.name)

    return orf_assign, orf_ambiguous, component_entry_idx, component_of


# ------------------------------------------ functional-consensus ancestral indel

def subtree_tip_sets(nodes_post):
    """{node.idx: frozenset(tip names under node)} in one postorder pass."""
    tips_under = {}
    for node in nodes_post:
        if node.is_tip:
            tips_under[node.idx] = frozenset((node.name,))
        else:
            acc = set()
            for child in node.children:
                acc |= tips_under[child.idx]
            tips_under[node.idx] = frozenset(acc)
    return tips_under


def mrca_node(carriers, root, tips_under):
    """Deepest node whose subtree contains every tip in `carriers`."""
    node = root
    while True:
        step = None
        for child in node.children:
            if carriers <= tips_under[child.idx]:
                step = child
                break
        if step is None:
            return node
        node = step


def functional_consensus_pin(states, functional_tips, root, tips_under, node_by_name, min_witnesses):
    """Return a {node_idx: PRESENT} pin map when an indel must be ancestral
    because it is shared by >= min_witnesses phylogenetically INDEPENDENT
    functional (complete-ORF) lineages; otherwise None.

    Rationale (real GUCY2F column-2772 case): complete-ORF lineages are
    trustworthy witnesses of the functional ancestral sequence. When several
    independent functional lineages all carry exactly the same indel,
    convergently gaining it in each is far less plausible than it being present
    in their common functional ancestor and then lost -- via random
    post-pseudogenization insertions/substitutions -- in the mostly-pseudogenized
    lineages that lack it.

    The indel is therefore fixed PRESENT not only at the functional carriers'
    common ancestor but at EVERY internal node on the paths from that ancestor
    down to each functional carrier. Pinning the MRCA alone is not enough: with
    the asymmetric cheap-gain cost the traceback would otherwise re-lose the gap
    just below the MRCA and re-gain it convergently on each functional branch --
    exactly the spurious independent deletion this rule exists to prevent. With
    the whole ancestor path pinned PRESENT, the functional lineages inherit the
    indel with no event and only the residue-carrier lineages record a loss.

    "Independent" means the functional carriers are not one clean clade: their
    MRCA subtree also contains a definite non-carrier (a tip scored ABSENT for
    this character). Two functional carriers that are each other's closest
    relatives (a lineage-specific shared indel in a single functional subclade)
    do NOT trigger the rule.
    """
    if not min_witnesses or min_witnesses < 2:
        return None
    carriers = frozenset(t for t in functional_tips if states.get(t) == PRESENT)
    if len(carriers) < min_witnesses:
        return None
    anc = mrca_node(carriers, root, tips_under)
    if not any(states.get(t) == ABSENT for t in tips_under[anc.idx]):
        return None
    pins = {}
    for tip_name in carriers:
        node = node_by_name[tip_name].parent
        while node is not None:
            pins[node.idx] = PRESENT
            if node is anc:
                break
            node = node.parent
    return pins or None


# ---------------------------------------------------------------- event merging

def merge_contiguous_same_type_indel_events(event_rows):
    """Real bug, found by inspecting real GUCA1C events directly: a single,
    contiguous, single-type indel on one branch (e.g. Nyctalus_aviator's real
    6bp deletion at columns 400-405) can be reported as two or more SEPARATE
    events (a 1bp deletion at 400, then a 5bp deletion at 401-405) whenever
    the breakpoint decomposition split that span into multiple characters --
    which it correctly does whenever some OTHER, unrelated lineage only
    shares PART of the span (here, Myotis_velifer independently shares just
    column 400), since each character is answering "is THIS exact span
    present/absent" as a comparable question across the whole tree. That
    per-character view is right for cross-species comparison, but for a
    single branch's own reported event list it fragments what is, from that
    branch's perspective alone, one real indel into an arbitrary-looking
    sequence of same-type fragments -- and can even manufacture spurious
    frameshift signal (e.g. a 1bp + 5bp fragment pair each individually not a
    multiple of 3, when their true combined 6bp span is in-frame).
    Directly requested: on ONE branch, two events of the SAME type (both
    "deletion", both "insertion_or_restoration", or both "ambiguous_indel_change")
    whose alignment spans are exactly contiguous (one starts the column
    immediately after the other ends, no gap) are the same real event and
    must be reported merged. Different types are never merged even when
    adjacent (a deletion right next to an insertion is two real, different
    events).

    Originally this excluded every shared_event=True row (an event whose
    origin is an internal/ancestral branch, inherited by multiple tips), on
    the reasoning that a shared event's origin branch is answered by
    ancestral parsimony rather than one lineage's own contiguous-span logic.
    Real data (CNGA3's Node95->Node94 branch: a confident 18bp deletion at
    466-483 immediately followed by a confident 3bp deletion at 484-486, both
    shared by the exact same six Myotis/Eptesicus tips) showed that reasoning
    was wrong -- the exact same breakpoint-decomposition fragmentation that
    splits a single-tip event also splits a shared, ancestral one, and it is
    just as fragmented and just as mergeable. Shared events ARE now merged,
    but only when every field that defines "the same reconstructed event" on
    that branch agrees, not just the branch and type: the exact same set of
    affected tips. That guards against the different failure mode visible on
    the same real branch, Node95->Node94, at columns 421 -- two GENUINELY
    DIFFERENT characters (one affecting 6 tips, the other a disjoint set of
    20) that happen to start at the same column and must stay two separate
    events, not merge just because they share a branch label.
    """
    MERGEABLE = {"deletion", "insertion_or_restoration", "ambiguous_indel_change"}
    clusterable, other = [], []
    for row in event_rows:
        if row["character_class"] == "indel" and row["biological_interpretation"] in MERGEABLE:
            clusterable.append(row)
        else:
            other.append(row)

    by_key = defaultdict(list)
    for row in clusterable:
        affected = frozenset(row["affected_tips"].split(","))
        by_key[(row["branch"], row["biological_interpretation"], affected)].append(row)

    merged = []
    for rows in by_key.values():
        # Real bug, found by inspecting real CNGA3 events directly: Myotis_
        # auriculus's own 466-483 and 484-486 deletions did NOT merge even
        # though they are genuinely contiguous, same branch, same type, same
        # affected tip -- because a THIRD, nested event on that same branch
        # (472-474, strictly inside 466-483, the same "independent interior
        # event" shape as GUCA1B's Nycteris case) sorted in BETWEEN them by
        # start column, so the naive "is the next start immediately after the
        # previous end" scan compared 466-483 against the nested 472-474
        # instead of against 484-486 and never saw the real adjacency at all.
        # A nested/contained interval must never participate in or interrupt
        # the sequential merge of the (non-nested) intervals that actually
        # partition the branch's real span -- it stays its own separate event,
        # exactly as decompose_run() already treats a genuine interior run.
        def is_nested(a, rows=rows):
            a_start, a_end = int(a["alignment_start"]), int(a["alignment_end"])
            for b in rows:
                if b is a:
                    continue
                b_start, b_end = int(b["alignment_start"]), int(b["alignment_end"])
                if b_start <= a_start and a_end <= b_end and (b_start, b_end) != (a_start, a_end):
                    return True
            return False
        nested = [r for r in rows if is_nested(r)]
        outer = [r for r in rows if not is_nested(r)]
        merged.extend(nested)
        outer.sort(key=lambda r: int(r["alignment_start"]))
        if outer:
            cluster = [outer[0]]
            for row in outer[1:]:
                if int(row["alignment_start"]) == int(cluster[-1]["alignment_end"]) + 1:
                    cluster.append(row)
                else:
                    merged.append(_merge_event_cluster(cluster))
                    cluster = [row]
            merged.append(_merge_event_cluster(cluster))
    return other + merged


def _merge_event_cluster(cluster):
    if len(cluster) == 1:
        return cluster[0]
    start = int(cluster[0]["alignment_start"])
    end = int(cluster[-1]["alignment_end"])
    length = end - start + 1
    widest = max(cluster, key=lambda r: int(r["event_length"]))  # most representative single character
    char_ids = [r["event_id"].split("|")[1] for r in cluster]
    kind_suffix = cluster[0]["event_id"].split("|")[2]
    merged = dict(cluster[0])
    merged.update({
        "event_id": f"{cluster[0]['gene']}|{'+'.join(char_ids)}|{kind_suffix}",
        "alignment_start": start,
        "alignment_end": end,
        "event_length": length,
        "length_mod_3": length % 3,
        "frame_effect": "in_frame" if length % 3 == 0 else "frameshift",
        "parsimony_score": format(sum(float(r["parsimony_score"]) for r in cluster), "g"),
        "delta_parsimony_support": min((r["delta_parsimony_support"] for r in cluster),
                                        key=lambda v: float("inf") if v in ("inf", "Infinity") else float(v)),
        "n_observed_present": widest["n_observed_present"],
        "n_observed_absent": widest["n_observed_absent"],
        "n_unknown": widest["n_unknown"],
        "observed_present_tips": widest["observed_present_tips"],
        "breakpoint_relationships": "merged_contiguous_same_type:" + "+".join(
            f"{r['event_id'].split('|')[1]}({r['alignment_start']}-{r['alignment_end']})" for r in cluster),
    })
    return merged


# ---------------------------------------------------------------- characters

def residue_span(sequence):
    """1-based (first, last) columns carrying a real ACGT residue, or None."""
    first = last = None
    for i, ch in enumerate(sequence, start=1):
        if ch in ACGT:
            if first is None:
                first = i
            last = i
    return (first, last) if first is not None else None


def event_is_terminal_incompleteness(occ, affected_tips, start, end):
    """True when an indel sits at the 5' or 3' end of its affected sequences.

    Such an indel reflects that the sequence is TRUNCATED (incomplete) at that
    end, not that the reading frame of the retained coding sequence is broken.
    Per the ORF-status rules a frameshift-length indel at the very beginning or
    end must NOT be treated as pseudogenizing evidence (real GRK7 case:
    Cynopterus_brachyotis carries only a 349 bp 3'-terminal deletion, no internal
    frameshift and no premature stop -- it is a partial sequence, not a
    pseudogene). An indel is terminal when NONE of its affected tips has any
    coding residue on one side of it (all residues lie to one side).
    """
    tips = [t for t in affected_tips if t in occ and occ[t] is not None]
    if not tips:
        return False
    has_5p = any(occ[t][0] <= start - 1 for t in tips)   # a residue strictly 5' of the indel
    has_3p = any(occ[t][1] >= end + 1 for t in tips)      # a residue strictly 3' of the indel
    return (not has_5p) or (not has_3p)


def gap_runs(sequence, gap="-"):
    """Maximal runs of gap in one sequence, as 1-based inclusive intervals."""
    runs, start = [], None
    for i, ch in enumerate(sequence, start=1):
        if ch == gap:
            if start is None:
                start = i
        elif start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(sequence)))
    return runs


def decompose_run(run, all_intervals, tolerance, _depth=0):
    """Split a gap run into a shared core plus its lineage-specific extension.

    Miniopterus australis carries 646-696 while natalensis and schreibersii
    carry 646-687.  Scored as raw runs these are three different characters and
    the shared 42 bp deletion disappears.  Scored as one 646-696 character the
    other two look like they gained 9 nt back.  Both are wrong.

    A run is therefore split against any OTHER observed run that it strictly
    contains AND with which it shares a breakpoint: 646-696 shares its left
    breakpoint with 646-687, so it becomes [646-687] + [688-696] - one shared
    deletion plus one extension, which is what the alignment actually shows.

    A contained run that shares NEITHER breakpoint is a separate, interior event
    (Nycteris thebaica 661-669 inside the Miniopterus deletion) and must not
    split anything.

    The shared core is decomposed RECURSIVELY, exactly like the residual, so a
    run is broken into the same atomic segments no matter which larger run it
    came from. Without this, decomposition is inconsistent: on real GRK7 data
    column 61 is an all-gap column, so most species carry a 61-66 gap while
    Molossus_molossus carries a longer 61-71 gap. (61,66) decomposes to
    (61,61)+(62,66) [it strictly contains ferrumequinum's 1 bp (61,61) run and
    shares that run's left breakpoint], but (61,71) used to pick its single
    largest one-breakpoint sub-run (61,66) as an undecomposed core, yielding
    (61,66)+(67,71). Molossus's (61,66) core then matched no other species'
    (61,66) segment -- because every other species had already been split into
    (61,61)+(62,66) -- so it clustered alone as a spurious singleton 61-66 indel
    character and was reported as a Molossus-only event, even though 61-66 is
    all-gap for Molossus and only 62-66 (ferrumequinum's insertion) carries any
    real residue. Recursively decomposing the core makes (61,71) resolve to
    (61,61)+(62,66)+(67,71), so Molossus joins the shared (62,66) character and
    only its genuine (67,71) extension is its own event. Interior runs that share
    neither breakpoint (Nycteris 661-669) still never split the core, so the
    documented Miniopterus/Nycteris behaviour is unchanged.
    """
    start, end = run
    if _depth > 64:
        return [(start, end)]
    best = None
    for cs, ce in all_intervals:
        if (cs, ce) == (start, end):
            continue
        if cs < start - tolerance or ce > end + tolerance:
            continue                                   # must be contained
        left = abs(cs - start) <= tolerance
        right = abs(ce - end) <= tolerance
        if left == right:                              # shares both, or neither
            continue
        if best is None or (ce - cs) > (best[1] - best[0]):
            best = (cs, ce)
    if best is None:
        return [(start, end)]
    cs, ce = best
    out = decompose_run((cs, ce), all_intervals, tolerance, _depth + 1)  # core, recursively
    if abs(cs - start) <= tolerance:                   # shared left -> residual right
        if ce + 1 <= end:
            out += decompose_run((ce + 1, end), all_intervals, tolerance, _depth + 1)
    else:                                              # shared right -> residual left
        if start <= cs - 1:
            out += decompose_run((start, cs - 1), all_intervals, tolerance, _depth + 1)
    return out


def cluster_indel_characters(alignment, tips, tolerance):
    """Group per-tip gap runs into breakpoint-defined characters.

    Two runs join the same character when BOTH breakpoints agree to within
    `tolerance` columns.  A run that shares only one breakpoint stays its own
    character but the relationship is recorded, which is what distinguishes a
    nested extension (Miniopterus australis 688-696 sharing a left breakpoint
    with the shared 646-687 deletion) from an independent indel.
    """
    raw = {tip: gap_runs(alignment[tip]) for tip in tips}
    all_intervals = sorted({r for runs in raw.values() for r in runs})
    segments_by_tip = {}
    for tip in tips:
        segs = []
        for run in raw[tip]:
            segs.extend(decompose_run(run, all_intervals, tolerance))
        segments_by_tip[tip] = sorted(set(segs))

    observed = [{"tip": tip, "start": s, "end": e}
                for tip in tips for s, e in segments_by_tip[tip]]

    clusters = []
    for run in sorted(observed, key=lambda r: (r["start"], r["end"], r["tip"])):
        for cluster in clusters:
            if (abs(cluster["start"] - run["start"]) <= tolerance
                    and abs(cluster["end"] - run["end"]) <= tolerance):
                cluster["tips"].add(run["tip"])
                cluster["members"].append(run)
                cluster["start"] = min(cluster["start"], run["start"])
                cluster["end"] = max(cluster["end"], run["end"])
                break
        else:
            clusters.append({"start": run["start"], "end": run["end"],
                             "tips": {run["tip"]}, "members": [run]})
    clusters.sort(key=lambda c: (c["start"], c["end"]))
    return clusters, segments_by_tip, raw


def breakpoint_relationships(clusters, tolerance):
    """Describe how each character relates to overlapping ones."""
    notes = defaultdict(list)
    for i, a in enumerate(clusters):
        for j, b in enumerate(clusters):
            if i == j:
                continue
            if b["end"] < a["start"] or b["start"] > a["end"]:
                continue
            same_start = abs(a["start"] - b["start"]) <= tolerance
            same_end = abs(a["end"] - b["end"]) <= tolerance
            if same_start and not same_end:
                kind = "shares_left_breakpoint_with"
            elif same_end and not same_start:
                kind = "shares_right_breakpoint_with"
            elif a["start"] <= b["start"] and a["end"] >= b["end"]:
                kind = "contains"
            elif a["start"] >= b["start"] and a["end"] <= b["end"]:
                kind = "nested_within"
            else:
                kind = "overlaps"
            notes[i].append(f"{kind}:IND{j + 1:04d}({b['start']}-{b['end']})")
    return notes


def indel_tip_states(alignment, tips, segments_by_tip, raw_runs, start, end, tolerance):
    """Score one breakpoint-defined indel character across the tips.

    This is complex indel coding, and the inapplicable case is the important
    one. A tip is PRESENT only when one of its OWN decomposed gap runs matches
    this character's breakpoints. If the entire character is gap in that tip
    only because a DIFFERENT larger deletion spans it, the state is UNKNOWN.
    If even one known residue occurs within the character, the full gap event is
    ABSENT -- so a smaller interior deletion does not make the taxon unknown for
    a larger character. On GUCA1B this makes Nycteris (661-669 deletion) ABSENT
    for the Miniopterus 646-687 event, while Miniopterus remain UNKNOWN for the
    smaller Nycteris character because their larger deletion spans all of it.
    """
    states = {}
    for tip in tips:
        matched = any(abs(s - start) <= tolerance and abs(e - end) <= tolerance
                      for s, e in segments_by_tip[tip])
        if matched:
            states[tip] = PRESENT
            continue
        seg = alignment[tip][start - 1:end]
        # A residue anywhere inside the character proves that this tip does NOT
        # carry the full breakpoint-defined gap event, even if it has a smaller
        # independent gap nested inside it. This is essential for GUCA1B:
        # Nycteris 661-669 is ABSENT for the larger Miniopterus 646-687 event.
        # Conversely, if the whole character is gap because a different larger
        # deletion spans it, the state is inapplicable/UNKNOWN rather than a
        # false carrier of the smaller event.
        if seg and any(c in ACGT for c in seg):
            states[tip] = ABSENT
        else:
            states[tip] = None
    return states


def stop_mask_characters(masked_stops, alignment, tips):
    """One character per exact mapped premature-STOP allele.

    Different STOP codons at the same aligned position are separate mutational
    characters. Raw STOPs whose MACSE frame phase is still shifted at that exact
    position are retained in diagnostics but are not reconstructed as independent
    nonsense mutations. A STOP downstream of compensated frameshifts (net MACSE
    frame correction mod 3 == 0) remains eligible as an independent event.
    """
    groups = defaultdict(set)
    alleles_at_span = defaultdict(dict)
    # Every tip's OWN raw-STOP evidence at each exact (species, start, end),
    # whether or not it passed validation, plus whether that occurrence was
    # SPECIFICALLY explained away as a frame-dependent consequence of an
    # earlier, unrelated frameshift (see the CNGA3/Phyllostomus discolor
    # regression). A tip with that specific, informative negative evidence
    # must never be resurrected as a carrier of another tip's validated
    # character merely because its alignment segment spells the same three
    # letters. A tip with no row here at all, or one that simply was never
    # independently classified for some other/unregistered reason, still
    # falls through to the ordinary segment-text-match check below -- see
    # test_shared_stop_not_missed_when_only_one_tip_is_registered, a real
    # PDE6H case where a genuinely shared mutation's second occurrence
    # legitimately never got its own independent classification at all.
    occurrence_by_tip_span = {}
    for row in masked_stops:
        start = safe_int(row.get("primary_alignment_start"))
        end = safe_int(row.get("primary_alignment_end"))
        species = row.get("species", "")
        codon = str(row.get("stop_codon", "")).upper()
        if start is None or end is None or not species or codon not in {"TAA", "TAG", "TGA"}:
            continue
        occurrence_by_tip_span[(species, start, end)] = {
            "codon": codon,
            "frame_shifted": trueish(row.get("frame_shifted_at_stop")),
        }
        if not trueish(row.get("pseudogenizing_event_candidate")):
            continue
        if "independent_stop_candidate" in row and not trueish(row.get("independent_stop_candidate")):
            continue
        groups[(start, end, codon)].add(species)
        alleles_at_span[(start, end)][species] = codon

    characters = []
    for (start, end, codon), carriers in sorted(groups.items()):
        states = {}
        for tip in tips:
            seg = alignment[tip][start - 1:end].upper()
            own = occurrence_by_tip_span.get((tip, start, end))
            if tip in carriers:
                states[tip] = PRESENT
            elif tip in alleles_at_span[(start, end)]:
                states[tip] = ABSENT  # another STOP allele is a different mutation
            elif own is not None and own["codon"] == codon and own["frame_shifted"]:
                # This tip's own raw STOP at this exact span/allele was
                # explicitly classified as a consequence of an earlier
                # frameshift, not an independent lesion. It must not be
                # resurrected as a carrier of this specific mutation just
                # because the alignment text matches.
                states[tip] = ABSENT
            elif seg and seg == codon:
                states[tip] = PRESENT
            elif seg and all(c in ACGT for c in seg):
                states[tip] = ABSENT
            else:
                states[tip] = None
        carriers_observed = set(carriers) | {t for t, s in states.items() if s == PRESENT}
        characters.append({"start": start, "end": end, "stop_codon": codon,
                           "tips": carriers_observed, "states": states})
    return characters


# ---------------------------------------------------------------- main

EVENT_HEADER = [
    "gene", "event_id", "character_class", "event_type", "biological_interpretation",
    "alignment_start", "alignment_end", "event_length", "length_mod_3", "frame_effect",
    "origin_node", "origin_is_tip", "parent_node", "branch",
    "shared_event", "n_affected_tips", "affected_tips",
    "reversal_below_origin", "secondary_changes_below_origin",
    "root_state", "parsimony_score", "delta_parsimony_support", "ambiguous_origin", "direction_confident",
    "parent_age", "child_age", "age_interval",
    "n_observed_present", "n_observed_absent", "n_unknown",
    "observed_present_tips", "terminal_incompleteness", "breakpoint_relationships", "coordinate_system",
]

CHAR_HEADER = [
    "gene", "character_id", "character_class", "alignment_start", "alignment_end",
    "character_length", "length_mod_3", "n_member_runs", "member_breakpoints",
    "n_observed_present", "n_observed_absent", "n_unknown", "root_state",
    "parsimony_score", "n_gains", "n_losses", "ambiguous_nodes", "history_ambiguous", "stop_codon",
    "functional_ancestral_indel", "breakpoint_relationships",
]

STATE_HEADER = ["gene", "character_id", "character_class", "node_label", "node_type", "state",
                "representative_state", "delta_parsimony_support", "ambiguous_at_node"]

FRAME_HEADER = ["gene", "node_label", "node_type", "net_indel_bp",
                "signed_net_length_change_from_root", "frame_shifted_by", "frame_offset_from_root",
                "frame_disrupted_provisional", "frame_currently_shifted",
                "frameshifting_events_in_history", "structural_state_ambiguous",
                "n_frameshift_indels_present", "n_stop_masks_present", "premature_stop_present",
                "frameshift_indel_ids", "stop_mask_ids"]


def main():
    parser = argparse.ArgumentParser(
        description="Breakpoint-coded indel/mask events placed on the dated tree by parsimony.")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--alignment", required=True,
                        help="NATIVE primary alignment (MACSE '!' rendered as '-')")
    parser.add_argument("--tree", required=True, help="Pruned dated tree")
    parser.add_argument("--masked-stops", default=None,
                        help="Step 01/02 masked in-frame STOP records (optional)")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dated", choices=["yes", "no"], default="yes")
    parser.add_argument("--gain-cost", type=float, default=0.35,
                        help="Cost of an ABSENT->PRESENT transition, for characters that can represent "
                             "a disabling lesion (a STOP allele, or a frameshifting indel) ONLY -- see the "
                             "orf_relevant check where this is used. Asymmetric with --loss-cost by "
                             "default: at a node with both intact and pseudogenized descendants, two "
                             "independent cheap gains (the gene broke separately in each lineage) is "
                             "treated as more plausible than one expensive loss (the ancestor was already "
                             "broken and one lineage reverted to functional) -- pseudogene-to-functional "
                             "reversal is biologically rare and should not tie 50:50 against a lesion "
                             "simply arising independently. Ordinary in-frame indels always use a neutral "
                             "1:1 cost regardless of this flag; set equal to --loss-cost to restore the "
                             "old symmetric 1:1 cost for disabling-lesion characters too.")
    parser.add_argument("--loss-cost", type=float, default=0.65,
                        help="Cost of a PRESENT->ABSENT transition for disabling-lesion characters "
                             "(a lesion reverting/apparently restored). Higher than --gain-cost by "
                             "default; see --gain-cost.")
    parser.add_argument("--tie-break", choices=["none", "ancestral", "terminal"], default="none",
                        help="Representative history for exact ties; ambiguity is always retained in output.")
    parser.add_argument("--breakpoint-tolerance", type=int, default=0,
                        help="Columns of slack when clustering indel breakpoints. "
                             "Start at 0 so real alignment disagreement is visible.")
    parser.add_argument("--min-carriers", type=int, default=2)
    parser.add_argument("--coordinate-system", default="primary_codon_alignment")
    parser.add_argument("--orf-status", default=None,
                        help="00_<gene>.orf_status.tsv (tip-level complete_orf calls). When given, enables "
                             "ORF-history-aware, pseudogenic-component-bounded character reconstruction "
                             "(Stage A/B/C): a provisional per-node Functional/Pseudogenic history is built "
                             "from this tip evidence ALONE (never from indel/STOP character calls, to avoid "
                             "circularity), and every confirmed F->P component's contribution to structural "
                             "character reconstruction ABOVE its own entry branch is capped at "
                             "--pseudogenic-boundary-cap-votes, so a densely-sampled pseudogenized clade "
                             "cannot out-vote a genuinely functional ancestor merely by tip count. Omit to "
                             "get the exact prior (v4.0 and earlier) uncapped behavior.")
    parser.add_argument("--orf-loss-cost", type=float, default=1.0,
                        help="Stage A: cost of an ordinary F->P (gene breaks) transition in the provisional "
                             "ORF-history reconstruction. Only the ratio to --orf-restoration-cost matters.")
    parser.add_argument("--orf-restoration-cost", type=float, default=2.0,
                        help="Stage A: cost of a P->F (apparently dead gene restored) transition -- kept "
                             "larger than --orf-loss-cost by default since true functional reversion "
                             "after pseudogenization is rare. A restoration:loss ratio sensitivity grid "
                             "(2/4/8/16) produced byte-identical real-data output (real GUCY2F, "
                             "Bat_genes_from_Song) and identical synthetic-test behavior at every ratio "
                             "tested; 2.0 (the smallest, most conservative ratio tested) is used as the "
                             "default per that stability -- see restoration_penalty_sensitivity.tsv.")
    parser.add_argument("--pseudogenic-boundary-cap-votes", type=float, default=2.0,
                        help="Stage C: maximum number of equivalent independent-tip votes a whole "
                             "pseudogenic component may contribute to its own parent's structural-character "
                             "cost, regardless of how many real descendant tips it contains -- the cap "
                             "itself is this value times max(gain_cost, loss_cost) for that character.")
    parser.add_argument("--min-functional-witnesses", type=int, default=2,
                        help="Functional-consensus ancestral-indel rule (requires --orf-status). If this "
                             "many or more phylogenetically INDEPENDENT complete-ORF (functional) lineages "
                             "carry exactly the same indel, that indel is fixed as PRESENT at their common "
                             "functional ancestor instead of being reconstructed as an independent gain on "
                             "each lineage -- convergent gain of the identical indel in several independent "
                             "functional lineages is biologically implausible; the parsimonious, biological "
                             "reading is an ancestral indel preserved in the functional lineages and lost "
                             "(via random post-pseudogenization mutations) in the pseudogenized ones (real "
                             "GUCY2F/column-2772 case). Set to 0 to disable; the minimum meaningful value "
                             "is 2.")
    args = parser.parse_args()

    if args.gain_cost <= 0 or args.loss_cost <= 0:
        raise SystemExit("--gain-cost and --loss-cost must be positive")
    if args.breakpoint_tolerance < 0:
        raise SystemExit("--breakpoint-tolerance must be >= 0")
    if args.orf_loss_cost <= 0 or args.orf_restoration_cost <= 0:
        raise SystemExit("--orf-loss-cost and --orf-restoration-cost must be positive")
    if args.pseudogenic_boundary_cap_votes <= 0:
        raise SystemExit("--pseudogenic-boundary-cap-votes must be positive")

    gene = args.gene
    dated = args.dated == "yes"
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    alignment, _order, aln_len = read_fasta(args.alignment)
    root = parse_newick(args.tree)
    root, collapsed = collapse_unifurcations(root)
    root, original_labels = apply_pensieve_labels(root)
    nodes_post = list(iter_nodes(root, "postorder"))
    tips = [n.name for n in nodes_post if n.is_tip]
    parent_of = {c.name: n.name for n in iter_nodes(root, "preorder") for c in n.children}

    missing = sorted(set(tips) - set(alignment))
    if missing:
        raise SystemExit(f"{len(missing)} tree tip(s) absent from the alignment: {missing[:10]}")
    extra = sorted(set(alignment) - set(tips))

    # ---- Stage A/B: provisional, tip-evidence-only ORF history and
    # pseudogenic-component boundaries (see pseudogenic_components()).
    # Computed once, before any structural character is reconstructed, and
    # reused for every character below -- never the other way around.
    component_entry_idx, component_of, orf_assign, orf_ambiguous = set(), {}, {}, set()
    functional_tips = frozenset()
    if args.orf_status:
        orf_status_rows = read_tsv(args.orf_status)
        tip_orf_states = read_tip_orf_states(orf_status_rows, tips)
        # complete-ORF (functional) tips: ABSENT in tip_orf_states means complete_orf.
        functional_tips = frozenset(t for t in tips if tip_orf_states.get(t) == ABSENT)
        orf_assign, orf_ambiguous, component_entry_idx, component_of = pseudogenic_components(
            nodes_post, root, tip_orf_states, args.orf_loss_cost, args.orf_restoration_cost, args.tie_break)
        print(f"[info] ORF-aware boundary evidence: {len(component_entry_idx)} pseudogenic component(s) "
              f"identified from tip-level ORF status; boundary cap = "
              f"{args.pseudogenic_boundary_cap_votes:g} equivalent tip vote(s)")
        write_tsv(
            [{"gene": gene, "node_label": n.name,
              "node_type": "tip" if n.is_tip else ("root" if n is root else "internal"),
              "provisional_orf_state": "ambiguous" if n.idx in orf_ambiguous else
                                       ("pseudogenic" if orf_assign.get(n.idx) == PRESENT else "functional"),
              "is_component_entry": n.idx in component_entry_idx,
              "pseudogenic_component": component_of.get(n.idx, "NA")}
             for n in nodes_post],
            out / f"03_{gene}.provisional_orf_history.tsv",
            ["gene", "node_label", "node_type", "provisional_orf_state", "is_component_entry",
             "pseudogenic_component"],
        )
        write_tsv(
            [{"gene": gene, "component_entry_node": n.name,
              "entry_branch": f"{n.parent.name}->{n.name}",
              "n_descendant_tips": sum(1 for t in iter_nodes(n, "preorder") if t.is_tip)}
             for n in iter_nodes(root, "preorder") if n.idx in component_entry_idx],
            out / f"03_{gene}.pseudogenic_components.tsv",
            ["gene", "component_entry_node", "entry_branch", "n_descendant_tips"],
        )

    ages = node_ages(root, dated)
    write_newick(root, out / f"03_{gene}.pensieve_labelled_dated_tree.nwk")
    label_map_rows = []
    for node in nodes_post:
        if node.is_tip:
            continue
        desc = sorted(t.name for t in iter_nodes(node, "preorder") if t.is_tip)
        label_map_rows.append({
            "gene": gene,
            "pensieve_node_label": node.name,
            "original_tree_label": original_labels.get(node.name, "NA"),
            "descendant_tip_count": len(desc),
            "descendant_tips": ",".join(desc),
            "rooted_descendant_tip_key": "|".join(desc),
            "node_age": ages.get(node.name, "NA"),
        })
    write_tsv(label_map_rows, out / f"03_{gene}.internode_label_map.tsv",
              ["gene", "pensieve_node_label", "original_tree_label",
               "descendant_tip_count", "descendant_tips",
               "rooted_descendant_tip_key", "node_age"])

    print(f"[info] {gene}: {len(tips)} tips, {len(nodes_post) - len(tips)} internal nodes, "
          f"alignment {aln_len} columns")
    if collapsed:
        print(f"[info] collapsed {len(collapsed)} degree-2 node(s)")
    if extra:
        print(f"[info] {len(extra)} alignment record(s) are not tree tips and were ignored")
    print(f"[info] parsimony gain={args.gain_cost:g} loss={args.loss_cost:g} "
          f"tie-break={args.tie_break} breakpoint-tolerance={args.breakpoint_tolerance}")

    # ---- build characters
    clusters, segments_by_tip, raw_runs = cluster_indel_characters(
        alignment, tips, args.breakpoint_tolerance)
    relationships = breakpoint_relationships(clusters, args.breakpoint_tolerance)
    characters = []
    for i, cluster in enumerate(clusters, start=1):
        characters.append({
            "id": f"IND{i:04d}",
            "klass": "indel",
            "start": cluster["start"],
            "end": cluster["end"],
            "states": indel_tip_states(alignment, tips, segments_by_tip, raw_runs,
                                       cluster["start"], cluster["end"],
                                       args.breakpoint_tolerance),
            "stop_codon": "NA",
            "members": cluster["members"],
            "relationships": ";".join(relationships.get(i - 1, [])) or "NA",
        })

    masked_stops = read_tsv(args.masked_stops) if args.masked_stops else []
    for i, character in enumerate(stop_mask_characters(masked_stops, alignment, tips), start=1):
        characters.append({
            "id": f"STOP{i:04d}",
            "klass": "stop_mask",
            "start": character["start"],
            "end": character["end"],
            "states": character["states"],
            "stop_codon": character.get("stop_codon", "NA"),
            "members": [],
            "relationships": "NA",
        })

    print(f"[info] characters: {sum(1 for c in characters if c['klass'] == 'indel')} indel, "
          f"{sum(1 for c in characters if c['klass'] == 'stop_mask')} stop-mask")

    # ---- reconstruct
    # Functional-consensus ancestral-indel rule (see functional_consensus_pin):
    # enabled only when tip ORF status is available and the witness threshold is
    # at least 2. tips_under is the subtree-tip map used to place an ancestral pin.
    use_functional_pin = bool(functional_tips) and args.min_functional_witnesses >= 2
    tips_under = subtree_tip_sets(nodes_post) if use_functional_pin else {}
    node_by_name = {n.name: n for n in nodes_post} if use_functional_pin else {}
    functional_ancestral_ids = []
    event_rows, char_rows, state_rows = [], [], []
    present_at = defaultdict(set)          # node -> set of character ids
    event_counter = 0

    for character in characters:
        states = character["states"]
        if not any(v == PRESENT for v in states.values()):
            continue
        length = character["end"] - character["start"] + 1
        # Asymmetric --gain-cost/--loss-cost is a biological prior about
        # PSEUDOGENIZATION specifically (functional->pseudogenized should not
        # tie 50:50 against pseudogenized->functional, the latter being much
        # rarer) -- it only means something for characters that can actually
        # represent a disabling lesion: a STOP allele, or a frameshifting
        # (length not a multiple of 3) indel. An ordinary in-frame indel's
        # gap-vs-residue reconstruction has no bearing on whether the ORF
        # stayed intact either way, so it always keeps the neutral 1:1 cost;
        # applying the disabling-lesion prior there too would just relabel
        # thousands of unrelated, harmless variants as spuriously "gap"-biased.
        orf_relevant = character["klass"] == "stop_mask" or (character["klass"] == "indel" and length % 3 != 0)
        gain_cost, loss_cost = (args.gain_cost, args.loss_cost) if orf_relevant else (1.0, 1.0)
        # Stage C: every structural character (not just ORF-relevant ones --
        # a pseudogenic clade's ordinary in-frame indels are just as
        # susceptible to the same many-descendants statistical inflation)
        # gets its pseudogenic-component contribution bounded, using THIS
        # character's own gain/loss costs to scale the cap.
        boundary_cap = (args.pseudogenic_boundary_cap_votes * max(gain_cost, loss_cost)
                        if component_entry_idx else None)
        # If two or more independent functional lineages share this indel, fix it
        # PRESENT at their common functional ancestor rather than let parsimony
        # reconstruct an implausible convergent gain on each functional branch.
        pinned = (functional_consensus_pin(states, functional_tips, root, tips_under,
                                           node_by_name, args.min_functional_witnesses)
                  if use_functional_pin and character["klass"] == "indel" else None)
        if pinned is not None:
            functional_ancestral_ids.append(character["id"])
        assign, score, ambiguous, delta = sankoff(
            nodes_post, states, gain_cost, loss_cost, args.tie_break,
            component_entry_idx=component_entry_idx, boundary_cap=boundary_cap, pinned=pinned)
        root_state = assign[root.idx]

        for node in nodes_post:
            is_ambiguous = node.idx in ambiguous
            if assign[node.idx] == PRESENT and not is_ambiguous:
                present_at[node.name].add(character["id"])
            if is_ambiguous:
                state_label = "ambiguous"
            elif character["klass"] == "indel":
                state_label = "gap" if assign[node.idx] == PRESENT else "residue"
            else:
                state_label = "stop_present" if assign[node.idx] == PRESENT else "stop_absent"
            state_rows.append({
                "gene": gene, "character_id": character["id"], "character_class": character["klass"],
                "node_label": node.name,
                "node_type": "tip" if node.is_tip else ("root" if node is root else "internal"),
                "state": state_label,
                "representative_state": "present" if assign[node.idx] == PRESENT else "absent",
                "delta_parsimony_support": format(delta.get(node.idx, 0.0), "g"),
                "ambiguous_at_node": is_ambiguous,
            })

        transitions = []
        for node in iter_nodes(root, "preorder"):
            if node.parent is None:
                continue
            if assign[node.parent.idx] == assign[node.idx]:
                continue
            transitions.append((node, "gain" if assign[node.idx] == PRESENT else "loss"))

        present_tips = {t for t, v in states.items() if v == PRESENT}
        absent_tips = {t for t, v in states.items() if v == ABSENT}
        unknown_tips = [t for t, v in states.items() if v is None]

        for node, kind in transitions:
            event_counter += 1
            subtree_tips = {t.name for t in iter_nodes(node, "preorder") if t.is_tip}
            # For a gain the affected tips carry the state; for a loss (a gap
            # lost = insertion/restoration) they are the ones that do NOT.
            affected = sorted((present_tips if kind == "gain" else absent_tips) & subtree_tips)
            secondary = [f"{n2.parent.name}->{n2.name}:{k2}"
                         for n2, k2 in transitions
                         if n2 is not node and n2.name in
                         {x.name for x in iter_nodes(node, "preorder")}]
            # Only ambiguity that can alter this branch transition/polarity
            # makes the event direction non-authoritative.  An unrelated
            # missing-data tie elsewhere in the tree is retained in the
            # character diagnostics but must not downgrade a clear origin.
            history_ambiguous = (
                root.idx in ambiguous or node.idx in ambiguous or node.parent.idx in ambiguous
            )
            direction_confident = not history_ambiguous
            if not direction_confident:
                interpretation = ("ambiguous_indel_change" if character["klass"] == "indel"
                                  else "ambiguous_stop_change")
            else:
                interpretation = (
                    ("deletion" if kind == "gain" else "insertion_or_restoration")
                    if character["klass"] == "indel"
                    else ("premature_stop_gained" if kind == "gain"
                          else "premature_stop_lost"))
            event_rows.append({
                "gene": gene,
                "event_id": f"{gene}|{character['id']}|{kind}{event_counter:04d}",
                "character_class": character["klass"],
                "event_type": f"{character['klass']}_{kind}",
                "biological_interpretation": interpretation,
                "alignment_start": character["start"],
                "alignment_end": character["end"],
                "event_length": length,
                "length_mod_3": length % 3,
                "frame_effect": ("in_frame" if length % 3 == 0 else "frameshift")
                                if character["klass"] == "indel" else "NA",
                "origin_node": node.name,
                "origin_is_tip": node.is_tip,
                "parent_node": node.parent.name,
                "branch": f"{node.parent.name}->{node.name}",
                "shared_event": (not node.is_tip) and len(affected) >= args.min_carriers,
                "n_affected_tips": len(affected),
                "affected_tips": ",".join(affected) or "NA",
                "reversal_below_origin": len(affected) < len(subtree_tips),
                "secondary_changes_below_origin": ";".join(secondary) or "NA",
                "root_state": "ambiguous" if root.idx in ambiguous else ("present" if root_state == PRESENT else "absent"),
                "parsimony_score": format(score, "g"),
                "delta_parsimony_support": format(delta.get(node.idx, 0.0), "g"),
                "ambiguous_origin": history_ambiguous or node.idx in ambiguous or node.parent.idx in ambiguous,
                "direction_confident": direction_confident,
                "parent_age": ages.get(node.parent.name, "NA"),
                "child_age": ages.get(node.name, "NA"),
                "age_interval": (f"{ages.get(node.parent.name)}-{ages.get(node.name)}"
                                 if dated else "NA"),
                "n_observed_present": len(present_tips),
                "n_observed_absent": len(absent_tips),
                "n_unknown": len(unknown_tips),
                "observed_present_tips": ",".join(sorted(present_tips)) or "NA",
                "breakpoint_relationships": character["relationships"],
                "coordinate_system": args.coordinate_system,
            })

        char_rows.append({
            "gene": gene, "character_id": character["id"],
            "character_class": character["klass"],
            "alignment_start": character["start"], "alignment_end": character["end"],
            "character_length": length, "length_mod_3": length % 3,
            "n_member_runs": len(character["members"]),
            "member_breakpoints": ";".join(
                sorted({f"{m['start']}-{m['end']}" for m in character["members"]})) or "NA",
            "n_observed_present": len(present_tips),
            "n_observed_absent": len(absent_tips),
            "n_unknown": len(unknown_tips),
            "root_state": "ambiguous" if root.idx in ambiguous else ("present" if root_state == PRESENT else "absent"),
            "parsimony_score": format(score, "g"),
            "n_gains": sum(1 for _, k in transitions if k == "gain"),
            "n_losses": sum(1 for _, k in transitions if k == "loss"),
            "ambiguous_nodes": ",".join(sorted(
                n.name for n in nodes_post if n.idx in ambiguous)) or "NA",
            "history_ambiguous": bool(ambiguous),
            "stop_codon": character.get("stop_codon", "NA"),
            "functional_ancestral_indel": pinned is not None,
            "breakpoint_relationships": character["relationships"],
        })

    event_rows = merge_contiguous_same_type_indel_events(event_rows)

    # Flag every event that is a 5'/3'-terminal indel (computed on the final,
    # post-merge span). A frameshift-length terminal indel means the sequence is
    # truncated there, not pseudogenized, so downstream ORF-status logic and the
    # figures exclude it from disabling evidence. stop_mask events never count.
    occ = {tip: residue_span(alignment[tip]) for tip in tips}
    for row in event_rows:
        if row.get("character_class") != "indel":
            row["terminal_incompleteness"] = False
            continue
        start = safe_int(row.get("alignment_start")); end = safe_int(row.get("alignment_end"))
        affected = [t for t in str(row.get("affected_tips", "")).split(",") if t and t != "NA"]
        row["terminal_incompleteness"] = bool(
            start is not None and end is not None
            and event_is_terminal_incompleteness(occ, affected, start, end))

    # ---- signed frame arithmetic: no ancestral nucleotide ASR required
    char_by_id = {c["id"]: c for c in characters}
    event_by_branch = defaultdict(list)
    for row in event_rows:
        event_by_branch[row["branch"]].append(row)

    frame_rows = []
    cumulative_delta = {root.name: 0}
    cumulative_frameshift_history = {root.name: 0}
    structural_ambiguous = {root.name: any(
        r["node_label"] == root.name and trueish(r["ambiguous_at_node"]) for r in state_rows
    )}
    for node in iter_nodes(root, "preorder"):
        if node.parent is not None:
            parent = node.parent.name
            branch = f"{parent}->{node.name}"
            delta_bp = 0
            fs_events = 0
            branch_ambiguous = False
            for ev in event_by_branch.get(branch, []):
                if ev["character_class"] != "indel":
                    continue
                if not trueish(ev.get("direction_confident")):
                    branch_ambiguous = True
                    continue
                length = int(ev["event_length"])
                if ev["biological_interpretation"] == "deletion":
                    delta_bp -= length
                elif ev["biological_interpretation"] == "insertion_or_restoration":
                    delta_bp += length
                if length % 3 != 0:
                    fs_events += 1
            cumulative_delta[node.name] = cumulative_delta[parent] + delta_bp
            cumulative_frameshift_history[node.name] = cumulative_frameshift_history[parent] + fs_events
            node_ambig = any(r["node_label"] == node.name and trueish(r["ambiguous_at_node"]) for r in state_rows)
            structural_ambiguous[node.name] = structural_ambiguous[parent] or branch_ambiguous or node_ambig

        ids = sorted(present_at.get(node.name, set()))
        indel_ids = [i for i in ids if char_by_id[i]["klass"] == "indel"]
        stop_ids = [i for i in ids if char_by_id[i]["klass"] == "stop_mask"]
        shifty = [i for i in indel_ids
                  if (char_by_id[i]["end"] - char_by_id[i]["start"] + 1) % 3]
        net = cumulative_delta.get(node.name, 0)
        offset = net % 3
        frame_rows.append({
            "gene": gene, "node_label": node.name,
            "node_type": "tip" if node.is_tip else ("root" if node is root else "internal"),
            "net_indel_bp": net,
            "signed_net_length_change_from_root": net,
            "frame_shifted_by": offset,
            "frame_offset_from_root": offset,
            "frame_disrupted_provisional": offset != 0,
            "frame_currently_shifted": offset != 0,
            "frameshifting_events_in_history": cumulative_frameshift_history.get(node.name, 0),
            "structural_state_ambiguous": structural_ambiguous.get(node.name, False),
            "n_frameshift_indels_present": len(shifty),
            "n_stop_masks_present": len(stop_ids),
            "premature_stop_present": bool(stop_ids),
            "frameshift_indel_ids": ",".join(shifty) or "NA",
            "stop_mask_ids": ",".join(stop_ids) or "NA",
        })

    write_tsv(event_rows, out / f"03_{gene}.alignment_events.tsv", EVENT_HEADER)
    write_tsv(char_rows, out / f"03_{gene}.alignment_characters.tsv", CHAR_HEADER)
    write_tsv(state_rows, out / f"03_{gene}.alignment_character_node_states.tsv", STATE_HEADER)
    write_tsv(frame_rows, out / f"03_{gene}.frame_arithmetic_by_node.tsv", FRAME_HEADER)

    shared = [r for r in event_rows if r["shared_event"]]
    ties = [r for r in event_rows if trueish(r["ambiguous_origin"])]
    print(f"[result] {len(event_rows)} events on {len({r['branch'] for r in event_rows})} branches; "
          f"{len(shared)} shared; {len(ties)} with an exact parsimony tie")
    if use_functional_pin:
        print(f"[info] functional-consensus rule (min {args.min_functional_witnesses} independent "
              f"complete-ORF witnesses): {len(functional_ancestral_ids)} indel character(s) fixed as "
              f"ancestral: {', '.join(functional_ancestral_ids) if functional_ancestral_ids else 'none'}")
    for row in sorted(shared, key=lambda r: (-int(r["event_length"]), int(r["alignment_start"])))[:15]:
        print(f"   {row['biological_interpretation']:<24} {row['alignment_start']}-{row['alignment_end']} "
              f"({row['event_length']} bp)  {row['branch']}  affects {row['n_affected_tips']}")
    print(f"[done] {out}")


if __name__ == "__main__":
    main()
