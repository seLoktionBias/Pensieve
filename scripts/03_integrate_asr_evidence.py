#!/usr/bin/env python3
"""Parse PAML marginal ASR and map it onto Pensieve internodes.

Pensieve v3.30 treats PAML's own ``rst`` declaration
``Nodes X to Y are ancestral`` as the source of truth for how many marginal
ancestral sequences that codeml run is expected to contain.  No n-2/n-1 guess is
used for ASR completeness.

The rooted user/Pensieve tree remains the biological authority.  With
``clock = 0`` the PAML TreeView Newick may contain one serialization/root vertex
in addition to the degree-3 biological unrooted internodes.  Its marginal
sequence is retained in a PAML-only audit FASTA but is never assigned to the
biological dated root.  Degree-3 PAML internodes are transferred to the rooted
reporting tree by exact root-independent tip tripartitions.

IndelMaP, when available, is concordance evidence only.  It never overwrites
Pensieve structural states or PAML nucleotide scaffolds.
"""
from __future__ import annotations

import argparse
import copy
import csv
import re
from collections import defaultdict
from pathlib import Path

from Bio import Phylo, SeqIO

DNA_GAPPED = set("ACGTN-")
UNKNOWN_PAML = set("?.*X!")


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
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=header, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path):
    records = {}
    order = []
    path = Path(path)
    if not path.exists():
        return records, order
    for record in SeqIO.parse(str(path), "fasta"):
        name = record.id.split()[0]
        records[name] = str(record.seq).upper().replace("U", "T")
        order.append(name)
    return records, order


def write_fasta(records, path, order=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if order is None:
        order = list(records)
    with path.open("w") as handle:
        for name in order:
            if name not in records:
                continue
            handle.write(f">{name}\n")
            seq = records[name]
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def safe_int(value):
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def trueish(value):
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def clean_name(name):
    return str(name or "").strip().strip("'\"")


def clean_paml_tip_name(name):
    return re.sub(r"^\d+_", "", clean_name(name), count=1)


def format_confidence(value):
    value = float(value)
    return str(int(value)) if value.is_integer() else format(value, "g")


def node_label_of(clade):
    if clade is None:
        return ""
    if clade.is_terminal():
        return clean_name(clade.name)
    label = clean_name(clade.name)
    if label:
        return label
    if clade.confidence is not None:
        return format_confidence(clade.confidence)
    return ""


def set_node_label(clade, label):
    clade.name = clean_name(label) or None
    clade.confidence = None


def write_newick(tree, path):
    safe_tree = copy.deepcopy(tree)
    for clade in safe_tree.get_nonterminals(order="preorder"):
        if clean_name(clade.name):
            clade.confidence = None
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(safe_tree, str(path), "newick", format_branch_length="%.17g")


def descendant_key(tips):
    return "|".join(sorted(clean_name(x) for x in tips))


def normalise_tree_tip_names(tree, paml_prefixes=False):
    names = []
    for tip in tree.get_terminals():
        tip.name = clean_paml_tip_name(tip.name) if paml_prefixes else clean_name(tip.name)
        names.append(tip.name)
    if any(not x for x in names):
        raise SystemExit("Tree contains an empty tip label")
    if len(names) != len(set(names)):
        duplicates = sorted(x for x in set(names) if names.count(x) > 1)
        raise SystemExit("Tree contains duplicated tip labels: " + ", ".join(duplicates[:20]))
    return tree


# ---------------------------------------------------------------------------
# PAML rst parsing.  These functions are used by BOTH the backend validator and
# downstream integration so there is only one interpretation of a real rst.


def parse_declared_ancestral_node_range(rst_file):
    """Return the inclusive ancestral-node range declared by PAML itself.

    We intentionally take the declaration immediately preceding the first
    marginal-ASR section, e.g. ``Nodes 104 to 205 are ancestral``.  This is the
    source of truth for ASR completeness; no expected count is derived from the
    number of tips.
    """
    lines = Path(rst_file).read_text(errors="replace").splitlines()
    marginal_i = next((i for i, line in enumerate(lines)
                       if "Marginal reconstruction of ancestral sequences" in line), None)
    if marginal_i is None:
        raise SystemExit("PAML rst lacks 'Marginal reconstruction of ancestral sequences'")

    matches = []
    pattern = re.compile(r"\bNodes?\s+(\d+)\s+to\s+(\d+)\s+are\s+ancestral\b", re.I)
    for i, line in enumerate(lines[:marginal_i + 1]):
        m = pattern.search(line)
        if m:
            matches.append((i, int(m.group(1)), int(m.group(2)), line.strip()))
    if not matches:
        raise SystemExit("PAML rst lacks the required 'Nodes X to Y are ancestral' declaration before marginal ASR")

    _i, start, end, raw = matches[-1]
    if end < start:
        raise SystemExit(f"Invalid PAML ancestral-node declaration: {raw}")
    numbers = list(range(start, end + 1))
    return {
        "start": start,
        "end": end,
        "count": len(numbers),
        "numbers": numbers,
        "labels": [f"PAML_Node{x}" for x in numbers],
        "raw_line": raw,
    }


def read_phylip_dimensions(path):
    with open(path) as handle:
        first = handle.readline().split()
    if len(first) < 2:
        raise SystemExit(f"Invalid PHYLIP header: {path}")
    return int(first[0]), int(first[1])


def normalise_paml_sequence(text):
    compact = re.sub(r"\s+", "", text.upper().replace("U", "T"))
    out = []
    for char in compact:
        if char in DNA_GAPPED:
            out.append(char)
        elif char in UNKNOWN_PAML:
            out.append("N")
    return "".join(out)


def parse_paml_marginal_sequences(rst_file, expected_length, expected_labels=None):
    """Parse only section (1), marginal reconstruction, from a PAML rst.

    If ``expected_labels`` is omitted it is taken directly from PAML's
    ``Nodes X to Y are ancestral`` declaration.  Every declared node must occur
    exactly once with an expected-length sequence.  Joint-ASR node records are
    deliberately ignored.
    """
    declaration = parse_declared_ancestral_node_range(rst_file)
    if expected_labels is None:
        expected_labels = declaration["labels"]
    expected_labels = list(expected_labels)
    expected_numbers = {int(x.replace("PAML_Node", "")) for x in expected_labels}

    lines = Path(rst_file).read_text(errors="replace").splitlines()
    marker_i = next(i for i, line in enumerate(lines)
                    if "Marginal reconstruction of ancestral sequences" in line)

    candidates_by_node = defaultdict(list)
    record_counts = defaultdict(int)
    current_node = None
    current_chunks = []

    def flush():
        nonlocal current_node, current_chunks
        if current_node is not None and current_chunks:
            seq = normalise_paml_sequence("".join(current_chunks))
            if seq:
                candidates_by_node[current_node].append(seq)
                record_counts[current_node] += 1
        current_node = None
        current_chunks = []

    for line in lines[marker_i + 1:]:
        stripped = line.strip()
        if "Joint reconstruction of ancestral sequences" in line:
            flush()
            break
        # PAML prints further material after the marginal DNA records: a
        # per-sequence/per-site accuracy summary ("Overall accuracy of the N
        # ancestral sequences:"), then "Amino acid sequences inferred by
        # codonml." with its own "Node #N  <one-letter AA string>" lines.
        # Those AA lines match the same "node #" pattern used for DNA records
        # below, so without an explicit stop here they get parsed as bogus
        # second/duplicate records for every node. None of this material, nor
        # anything from a repeated tree-view block, is sequence continuation.
        if any(token in line for token in [
            "tree with node labels for Rod Page's TreeView",
            "Overall accuracy of the",
            "Amino acid sequences inferred by codonml",
        ]):
            flush()
            break
        match = re.match(r"^\s*node\s*#?\s*(\d+)\s+(.+?)\s*$", line, re.I)
        if match:
            flush()
            current_node = int(match.group(1))
            current_chunks = [match.group(2)]
            continue
        if current_node is not None and stripped and re.fullmatch(r"[A-Za-z?*!.\-\s]+", line):
            current_chunks.append(line)
            continue
        if current_node is not None and not stripped:
            flush()
    flush()

    result = {}
    diagnostics = []
    duplicate_nodes = []
    for node_number in sorted(candidates_by_node):
        sequences = candidates_by_node[node_number]
        unique = []
        for seq in sequences:
            if seq not in unique:
                unique.append(seq)
        exact = [seq for seq in unique if len(seq) == expected_length and set(seq).issubset(DNA_GAPPED)]
        label = f"PAML_Node{node_number}"
        if len(sequences) > 1:
            duplicate_nodes.append(label)
        if exact:
            result[label] = exact[0]
        diagnostics.append({
            "paml_node_label": label,
            "records_in_marginal_section": len(sequences),
            "candidate_lengths": ",".join(str(len(x)) for x in unique) if unique else "NA",
            "selected_length": len(result.get(label, "")),
        })

    observed_numbers = {int(x.replace("PAML_Node", "")) for x in result}
    missing_numbers = sorted(expected_numbers - observed_numbers)
    extra_numbers = sorted(observed_numbers - expected_numbers)
    if duplicate_nodes:
        raise SystemExit("Duplicate marginal node records in PAML rst: " + ", ".join(duplicate_nodes[:30]))
    if missing_numbers:
        raise SystemExit(
            "Missing/invalid expected-length marginal sequences for PAML nodes: "
            + ", ".join(f"PAML_Node{x}" for x in missing_numbers[:30])
        )
    if extra_numbers:
        raise SystemExit(
            "Marginal ASR contains node(s) outside PAML's declared ancestral range: "
            + ", ".join(f"PAML_Node{x}" for x in extra_numbers[:30])
        )
    if len(result) != declaration["count"]:
        raise SystemExit(
            f"Marginal ASR count mismatch against PAML declaration: observed={len(result)} "
            f"declared={declaration['count']}"
        )
    return result, diagnostics, declaration


def validate_paml_marginal_output(rst_file, phylip_file):
    """Validate exactly the marginal product Pensieve consumes downstream."""
    rst_file = Path(rst_file)
    phylip_file = Path(phylip_file)
    if not rst_file.exists() or rst_file.stat().st_size == 0:
        raise SystemExit(f"missing_or_empty_rst:{rst_file}")
    if not phylip_file.exists() or phylip_file.stat().st_size == 0:
        raise SystemExit(f"missing_or_empty_phylip:{phylip_file}")
    tip_count, alignment_length = read_phylip_dimensions(phylip_file)
    marginal, diagnostics, declaration = parse_paml_marginal_sequences(rst_file, alignment_length)
    text = rst_file.read_text(errors="replace")
    row = {
        "marginal_asr_valid": True,
        "reason": "all_rst_declared_marginal_nodes_present_once_at_expected_alignment_length",
        "tip_count": tip_count,
        "declared_ancestral_start": declaration["start"],
        "declared_ancestral_end": declaration["end"],
        "expected_internal_sequences": declaration["count"],
        "observed_internal_sequences": len(marginal),
        "alignment_length": alignment_length,
        "joint_section_present": "Joint reconstruction of ancestral sequences" in text,
        "overall_accuracy_marker_present": "Overall accuracy of the" in text,
        "validation_basis": "rst_line:Nodes_X_to_Y_are_ancestral",
    }
    return row, marginal, diagnostics, declaration


def reconstruct_missing_paml_serialization_vertex_label(tree, declared_labels):
    """Label an omitted PAML Newick serialization root from the rst declaration.

    PAML sometimes omits the label on the outer Newick/root vertex while still
    declaring/reconstructing that node.  If exactly one declared ancestral label
    is absent from the parsed internal labels, assign that single missing label
    to the Newick root.  No n-derived numbering guess is made.
    """
    existing = node_label_of(tree.root)
    if existing:
        return existing, "paml_treeview_label"
    declared = set(declared_labels)
    observed = {node_label_of(c) for c in tree.get_nonterminals(order="preorder") if node_label_of(c)}
    missing = sorted(declared - observed)
    if len(missing) == 1:
        set_node_label(tree.root, missing[0])
        return missing[0], "reconstructed_from_unique_rst_declared_label_missing_at_newick_root"
    raise SystemExit(
        "PAML TreeView root is unlabelled and cannot be assigned uniquely from the rst ancestral-node declaration: "
        f"missing_declared_labels={','.join(missing[:20]) or 'none'}"
    )


def extract_paml_labelled_tree(rst_file, declared_labels=None):
    text = Path(rst_file).read_text(errors="replace")
    marker = "tree with node labels for Rod Page's TreeView"
    pos = text.find(marker)
    if pos < 0:
        raise SystemExit(f"PAML rst lacks required marker: {marker}")
    start = text.find("(", pos)
    end = text.find(";", start)
    if start < 0 or end < 0:
        raise SystemExit("Could not extract the PAML labelled Newick tree after the TreeView marker")
    raw = text[start:end + 1]
    raw = re.sub(r"(?<=[(,])\s*\d+_", "", raw)
    raw = re.sub(r"\)\s*(\d+)\s*(?=[:,);])", r")PAML_Node\1", raw)
    raw = re.sub(r"\s+", "", raw)
    tmp = Path(rst_file).with_name(".pensieve_paml_labelled_tree.tmp.nwk")
    tmp.write_text(raw + "\n")
    try:
        tree = Phylo.read(str(tmp), "newick")
    finally:
        tmp.unlink(missing_ok=True)
    normalise_tree_tip_names(tree, paml_prefixes=False)
    if declared_labels is None:
        declared_labels = parse_declared_ancestral_node_range(rst_file)["labels"]
    serialization_label, serialization_source = reconstruct_missing_paml_serialization_vertex_label(tree, declared_labels)
    unnamed = [c for c in tree.get_nonterminals() if not node_label_of(c)]
    if unnamed:
        raise SystemExit(f"PAML labelled tree contains {len(unnamed)} unnamed internal nodes")
    labels = [node_label_of(c) for c in tree.get_nonterminals()]
    if len(labels) != len(set(labels)):
        raise SystemExit("PAML labelled tree contains duplicated internal labels")
    undeclared = sorted(set(labels) - set(declared_labels))
    if undeclared:
        raise SystemExit("PAML TreeView contains internal labels outside the declared ancestral range: " + ", ".join(undeclared[:20]))
    return tree, raw, serialization_label, serialization_source


# ---------------------------------------------------------------------------
# Tree topology and mapping helpers.


def undirected_tree_graph(tree):
    adjacency = defaultdict(list)
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            adjacency[parent].append((child, child.branch_length))
            adjacency[child].append((parent, child.branch_length))
    return adjacency


def component_tip_set(adjacency, start, blocked):
    tips = set()
    stack = [(start, blocked)]
    while stack:
        node, previous = stack.pop()
        if node.is_terminal():
            tips.add(clean_name(node.name))
        for neighbour, _length in adjacency[node]:
            if neighbour is previous:
                continue
            stack.append((neighbour, node))
    return tips


def node_partition_signature(tree, clade, adjacency=None):
    if adjacency is None:
        adjacency = undirected_tree_graph(tree)
    partitions = [
        descendant_key(component_tip_set(adjacency, neighbour, clade))
        for neighbour, _length in adjacency[clade]
    ]
    return "||".join(sorted(partitions))


def internal_signature_index(tree, require_degree_three=True):
    adjacency = undirected_tree_graph(tree)
    index = defaultdict(list)
    for clade in tree.get_nonterminals(order="preorder"):
        if require_degree_three and len(adjacency[clade]) != 3:
            continue
        index[node_partition_signature(tree, clade, adjacency)].append(clade)
    return index, adjacency


def rooted_internal_index(tree, source_name, require_labels=False):
    index = {}
    labels = {}
    for clade in tree.get_nonterminals(order="postorder"):
        key = descendant_key(t.name for t in clade.get_terminals())
        if key in index:
            raise SystemExit(f"{source_name} contains duplicate rooted internal descendant sets: {key}")
        label = node_label_of(clade)
        if require_labels and not label:
            raise SystemExit(f"{source_name} contains an unlabelled internal node for descendant set: {key}")
        index[key] = clade
        labels[key] = label
    return index, labels


def topology_metrics(tree):
    graph = undirected_tree_graph(tree)
    signatures, _ = internal_signature_index(tree, require_degree_three=True)
    return {
        "tip_set": {clean_name(t.name) for t in tree.get_terminals()},
        "signature_set": set(signatures),
        "root_degree": len(graph[tree.root]),
    }


def validate_same_unrooted_topology(reference, target, target_name):
    ref = topology_metrics(reference)
    obs = topology_metrics(target)
    if ref["tip_set"] != obs["tip_set"]:
        raise SystemExit(
            f"Tip-set mismatch for {target_name}: only_user={','.join(sorted(ref['tip_set']-obs['tip_set'])[:20]) or 'NA'}; "
            f"only_{target_name}={','.join(sorted(obs['tip_set']-ref['tip_set'])[:20]) or 'NA'}"
        )
    if ref["signature_set"] != obs["signature_set"]:
        raise SystemExit(f"Unrooted topology mismatch between user tree and {target_name}")


def normalise_tripartition_signature(value):
    value = str(value or "").strip()
    return value[len("TRIPARTITION::"):] if value.startswith("TRIPARTITION::") else value


def unique_internal_tripartitions(tree, source_name, require_labels=False):
    """Index degree-3 biological internodes by root-independent tip tripartition."""
    raw_index, graph = internal_signature_index(tree, require_degree_three=True)
    index = {}
    labels = {}
    for raw_signature, clades in raw_index.items():
        signature = normalise_tripartition_signature(raw_signature)
        if len(clades) != 1:
            raise SystemExit(f"{source_name} has a non-unique internal tripartition: {signature}")
        clade = clades[0]
        label = node_label_of(clade)
        if require_labels and not label:
            raise SystemExit(f"{source_name} has an unlabelled biological internode: {signature}")
        index[signature] = clade
        labels[signature] = label
    expected = len(tree.get_terminals()) - 2
    if len(index) != expected:
        degrees = [len(graph[c]) for c in tree.get_nonterminals(order="preorder")]
        raise SystemExit(
            f"{source_name} has {len(index)} degree-3 biological internodes; expected {expected}. "
            f"Internal degrees={','.join(map(str, degrees[:30]))}"
        )
    return index, labels


def transfer_paml_labels_to_dated_tree_by_tripartition(paml_tree, reporting_tree):
    reporting = copy.deepcopy(reporting_tree)
    user_index, _ = unique_internal_tripartitions(reporting, "dated reporting tree")
    paml_index, paml_labels = unique_internal_tripartitions(paml_tree, "PAML TreeView tree", require_labels=True)
    if set(user_index) != set(paml_index):
        raise SystemExit("PAML internodes cannot be transferred to the dated tree: tripartitions differ")
    root_original = node_label_of(reporting.root) or "NA"
    set_node_label(reporting.root, "UserRoot")
    rows = [{
        "paml_node_label": "NA",
        "reporting_tree_original_node_label": root_original,
        "reporting_tree_node_label": "UserRoot",
        "is_user_tree_root": True,
        "descendant_tip_count": len(reporting.get_terminals()),
        "descendant_tips": ",".join(sorted(t.name for t in reporting.get_terminals())),
        "rooted_descendant_tip_key": descendant_key(t.name for t in reporting.get_terminals()),
        "root_independent_node_signature": "ROOT_EDGE",
        "mapping_method": "biological_root_not_assigned_clock0_paml_serialization_sequence",
        "mapping_status": "NO_PAML_LABEL_EXPECTED",
    }]
    for signature, target in user_index.items():
        original = node_label_of(target) or "NA"
        paml_label = paml_labels[signature]
        set_node_label(target, paml_label)
        rows.append({
            "paml_node_label": paml_label,
            "reporting_tree_original_node_label": original,
            "reporting_tree_node_label": paml_label,
            "is_user_tree_root": False,
            "descendant_tip_count": len(target.get_terminals()),
            "descendant_tips": ",".join(sorted(t.name for t in target.get_terminals())),
            "rooted_descendant_tip_key": descendant_key(t.name for t in target.get_terminals()),
            "root_independent_node_signature": "TRIPARTITION::" + signature,
            "mapping_method": "exact_root_independent_three_tip_partition",
            "mapping_status": "PASS",
        })
    return reporting, rows


# ---------------------------------------------------------------------------
# Optional IndelMaP concordance.


def indelmap_concordance(gene, r2, out, pensieve_tree):
    # Matches indelMaP_ASR.py's own naming: it appends fixed suffixes to
    # --output_file (e.g. "<output_file>_tree.nwk",
    # "<output_file>_internal_ancestral_reconstruction.fas"); it does not
    # insert an "_ASR" infix. Verified against a real run of external/indelMaP.
    tree_file = r2 / "indelmap_asr" / f"{gene}.indelmap_tree.nwk"
    fasta_file = r2 / "indelmap_asr" / f"{gene}.indelmap_internal_ancestral_reconstruction.fas"
    status_file = out / f"03_{gene}.indelmap_concordance_status.tsv"
    concord_file = out / f"03_{gene}.indelmap_concordance.tsv"
    header = ["gene", "character_id", "pensieve_node_label", "indelmap_node_label", "pensieve_state", "indelmap_state", "agreement"]
    if not tree_file.exists() or not fasta_file.exists() or tree_file.stat().st_size == 0 or fasta_file.stat().st_size == 0:
        write_tsv([{"gene": gene, "status": "not_available", "details": "IndelMaP disabled, failed, or did not emit expected ancestral files"}], status_file)
        write_tsv([], concord_file, header)
        return

    indel_tree = normalise_tree_tip_names(Phylo.read(str(tree_file), "newick"))
    validate_same_unrooted_topology(pensieve_tree, indel_tree, "IndelMaP tree")
    pens_index, pens_labels = rooted_internal_index(pensieve_tree, "Pensieve tree", require_labels=True)
    indel_index, indel_labels = rooted_internal_index(indel_tree, "IndelMaP tree", require_labels=True)
    if set(pens_index) != set(indel_index):
        write_tsv([{"gene": gene, "status": "rooted_topology_disagreement", "details": "IndelMaP rooted descendant sets differ from Pensieve tree"}], status_file)
        write_tsv([], concord_file, header)
        return

    indel_sequences, _ = read_fasta(fasta_file)
    characters = {r["character_id"]: r for r in read_tsv(out / f"03_{gene}.alignment_characters.tsv") if r.get("character_class") == "indel"}
    pstate = {(r["character_id"], r["node_label"]): r.get("state", "unknown") for r in read_tsv(out / f"03_{gene}.alignment_character_node_states.tsv")}
    rows = []
    for key in sorted(pens_index):
        plabel = pens_labels[key]
        ilabel = indel_labels[key]
        seq = indel_sequences.get(ilabel)
        if seq is None:
            continue
        for cid, char in characters.items():
            start = safe_int(char.get("alignment_start")); end = safe_int(char.get("alignment_end"))
            if start is None or end is None or end > len(seq):
                istate = "unknown"
            else:
                segment = seq[start - 1:end].upper()
                if segment and all(c == "-" for c in segment):
                    istate = "gap"
                elif segment and all(c in "ACGTN" for c in segment) and any(c in "ACGT" for c in segment):
                    istate = "residue"
                else:
                    istate = "unknown"
            ps = pstate.get((cid, plabel), "unknown")
            comparable = ps in {"gap", "residue"} and istate in {"gap", "residue"}
            rows.append({
                "gene": gene, "character_id": cid, "pensieve_node_label": plabel,
                "indelmap_node_label": ilabel, "pensieve_state": ps, "indelmap_state": istate,
                "agreement": (ps == istate) if comparable else "NA",
            })
    write_tsv(rows, concord_file, header)
    comparable = [r for r in rows if r["agreement"] != "NA"]
    disagree = [r for r in comparable if not trueish(r["agreement"])]
    write_tsv([{
        "gene": gene, "status": "available", "n_comparable_states": len(comparable),
        "n_disagreements": len(disagree), "details": "IndelMaP is concordance evidence only; it never overwrites Pensieve states",
    }], status_file)


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Parse codeml marginal ASR and map PAML internodes to Pensieve nodes.")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--results02-dir", required=True)
    parser.add_argument("--results00-dir", required=True)  # retained for CLI compatibility; not used
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dated", choices=["yes", "no"], default="yes")
    parser.add_argument("--on-missing-root-sequence", choices=["fail", "warn"], default="warn")
    args = parser.parse_args()

    gene = args.gene
    r2 = Path(args.results02_dir)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    rst = r2 / "paml_codon_asr" / "rst"
    phylip = r2 / f"02_{gene}.codon_for_paml.phy"
    validation, paml_all, diagnostics, declaration = validate_paml_marginal_output(rst, phylip)
    alignment_length = validation["alignment_length"]
    write_tsv(diagnostics, out / f"03_{gene}.paml_marginal_asr_parse_diagnostics.tsv")

    paml_tree, _raw, serialization_label, serialization_source = extract_paml_labelled_tree(
        rst, declaration["labels"]
    )
    write_newick(paml_tree, out / f"03_{gene}.paml_tree_from_rst.nwk")

    reporting_tree = normalise_tree_tip_names(Phylo.read(str(r2 / f"02_{gene}.tree_for_asr.nwk"), "newick"))
    validate_same_unrooted_topology(reporting_tree, paml_tree, "PAML TreeView tree")

    # PAML may reconstruct one extra sequence for the degree-2 Newick
    # serialization/root vertex.  It is a valid PAML marginal sequence and is
    # retained for audit, but it is not the biological dated root sequence.
    _, biological_paml_labels = unique_internal_tripartitions(paml_tree, "PAML TreeView tree", require_labels=True)
    biological_label_set = set(biological_paml_labels.values())
    all_label_set = set(paml_all)
    extra_paml_labels = sorted(all_label_set - biological_label_set)
    missing_biological = sorted(biological_label_set - all_label_set)
    if missing_biological:
        raise SystemExit("PAML marginal ASR lacks biological internode sequence(s): " + ", ".join(missing_biological[:20]))
    if len(extra_paml_labels) > 1:
        raise SystemExit(
            "More than one rst-declared PAML ancestral sequence does not map to a degree-3 biological internode: "
            + ", ".join(extra_paml_labels[:20])
        )
    if extra_paml_labels and extra_paml_labels[0] != serialization_label:
        raise SystemExit(
            "The single non-biological PAML marginal sequence is not the TreeView serialization/root label: "
            f"extra={extra_paml_labels[0]} serialization={serialization_label}"
        )
    biological_marginal = {label: paml_all[label] for label in biological_label_set}

    # Audit all PAML-declared marginal sequences separately from the biological
    # non-root subset used by Pensieve.
    numeric_order = [f"PAML_Node{x}" for x in declaration["numbers"]]
    write_fasta(paml_all, out / f"03_{gene}.paml_marginal_asr_all_declared_nodes.fa", numeric_order)

    labelled_reporting, mapping_rows = transfer_paml_labels_to_dated_tree_by_tripartition(paml_tree, reporting_tree)
    paml_phylo_order = [
        node_label_of(c) for c in labelled_reporting.find_clades(order="preorder")
        if (not c.is_terminal()) and c is not labelled_reporting.root
    ]
    if set(paml_phylo_order) != set(biological_marginal):
        raise SystemExit("Cannot phylogenetically order PAML biological marginal ASR: mapped node-label sets differ")
    write_fasta(biological_marginal, out / f"03_{gene}.paml_marginal_asr.fa", paml_phylo_order)
    write_newick(labelled_reporting, out / f"03_{gene}.paml_labeled_reporting_tree.nwk")
    write_tsv(mapping_rows, out / f"03_{gene}.paml_to_reporting_tree_node_map.tsv")

    pensieve_file = out / f"03_{gene}.pensieve_labelled_dated_tree.nwk"
    if not pensieve_file.exists():
        raise SystemExit(f"Missing Pensieve event tree; run the events stage first: {pensieve_file}")
    pensieve_tree = normalise_tree_tip_names(Phylo.read(str(pensieve_file), "newick"))
    validate_same_unrooted_topology(reporting_tree, pensieve_tree, "Pensieve event tree")

    pidx, plabels = rooted_internal_index(pensieve_tree, "Pensieve tree", require_labels=True)
    _ridx, rlabels = rooted_internal_index(labelled_reporting, "PAML-labelled reporting tree", require_labels=True)
    root_key = descendant_key(t.name for t in pensieve_tree.root.get_terminals())
    crosswalk = []
    for key in sorted(pidx, key=lambda k: (len(k.split("|")), k)):
        pens_label = plabels[key]
        paml_label = "NA" if key == root_key else rlabels.get(key, "NA")
        crosswalk.append({
            "gene": gene,
            "pensieve_node_label": pens_label,
            "paml_node_label": paml_label,
            "is_user_tree_root": key == root_key,
            "rooted_descendant_tip_key": key,
            "mapping_basis": (
                "biological_root_not_assigned_clock0_paml_serialization_sequence"
                if key == root_key else
                "exact_rooted_descendant_tip_set_after_tripartition_transfer"
            ),
            "mapping_status": "NO_PAML_SEQUENCE_EXPECTED" if key == root_key else ("PASS" if paml_label != "NA" else "FAIL"),
        })
    if any(r["mapping_status"] == "FAIL" for r in crosswalk):
        raise SystemExit("At least one non-root Pensieve internode could not be mapped to a PAML marginal sequence")
    write_tsv(crosswalk, out / f"03_{gene}.internode_label_crosswalk.tsv", [
        "gene", "pensieve_node_label", "paml_node_label", "is_user_tree_root",
        "rooted_descendant_tip_key", "mapping_basis", "mapping_status",
    ])

    write_tsv([{
        "gene": gene,
        "tip_count": len(reporting_tree.get_terminals()),
        "rst_declared_ancestral_start": declaration["start"],
        "rst_declared_ancestral_end": declaration["end"],
        "rst_declared_marginal_sequence_count": declaration["count"],
        "parsed_marginal_sequence_count": len(paml_all),
        "mapped_biological_nonroot_sequence_count": len(biological_marginal),
        "serialization_only_sequence_count": len(extra_paml_labels),
        "serialization_only_paml_label": extra_paml_labels[0] if extra_paml_labels else "NA",
        "alignment_length": alignment_length,
        "paml_newick_serialization_vertex_label": serialization_label,
        "paml_newick_serialization_vertex_label_source": serialization_source,
        "unrooted_topology_matches_user_tree": True,
        "biological_root_sequence_available_from_clock0_paml": False,
        "asr_validation_basis": "rst_line:Nodes_X_to_Y_are_ancestral",
        "status": "PASS",
    }], out / f"03_{gene}.tree_topology_and_node_count_audit.tsv")

    indelmap_concordance(gene, r2, out, pensieve_tree)
    warning_lines = [
        "PAML marginal-ASR completeness was validated from its own 'Nodes X to Y are ancestral' declaration.",
        "The biological dated root is not assigned the clock=0 PAML Newick serialization/root sequence; root ORF status remains unavailable.",
    ]
    if extra_paml_labels:
        warning_lines.append(
            f"PAML marginal sequence {extra_paml_labels[0]} belongs to the TreeView serialization/root vertex and is retained only in the all-declared-nodes audit FASTA."
        )
    (out / f"03_{gene}.warnings.txt").write_text("\n".join(warning_lines) + "\n")
    print(
        f"Finished PAML integration for {gene}; rst_declared_marginal={len(paml_all)}; "
        f"biological_nonroot_mapped={len(biological_marginal)}; alignment_length={alignment_length}"
    )


if __name__ == "__main__":
    main()
