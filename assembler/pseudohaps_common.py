#!/usr/bin/env python3
"""
pseudohaps_common.py

Shared GFA graph engine and assembled-chain discovery reused by the anchor-based
pseudohaplotype scripts (pseudohaps_siblings.py, pseudohaps_build.py,
pseudohaps_anchors.py, pseudohaps_regional.py). Also provides a legacy
standalone graph-topology pseudohaplotype pipeline via main().

Analyze GFA files and identify assembled chains by finding paths
from start nodes (no incoming edges) to end nodes (no outgoing edges).
"""

from collections import defaultdict
import os
import time
from typing import Dict, List, Set, Tuple, Optional


class GFAGraph:
    """Data structure to represent a GFA graph."""

    def __init__(self):
        self.nodes: Dict[str, str] = {}
        self.out_edges: Dict[str, List[Tuple[str, str, str, int]]] = defaultdict(list)
        self.in_edges: Dict[str, List[Tuple[str, str, str, int]]] = defaultdict(list)
        self.all_nodes: Set[str] = set()
        # Store edge-specific overlaps for path sequence construction
        self.edge_overlaps: Dict[Tuple[str, str], int] = {}

    def add_node(self, node_id: str, sequence: str):
        self.nodes[node_id] = sequence
        self.all_nodes.add(node_id)

    def add_edge(self, from_node: str, from_orient: str, to_node: str, to_orient: str, overlap: int = 0):
        self.out_edges[from_node].append((to_node, from_orient, to_orient, overlap))
        self.in_edges[to_node].append((from_node, from_orient, to_orient, overlap))
        self.all_nodes.add(from_node)
        self.all_nodes.add(to_node)
        # Store edge-specific overlap
        self.edge_overlaps[(from_node, to_node)] = overlap
    
    def resolve_overlaps(self) -> None:
        """
        Resolve overlaps by trimming node sequences.
        
        Assumption: All outgoing edges from a parent have the same overlap.
        
        For each node:
        1. If it has outgoing edges: trim suffix by the overlap value
           (all outgoing edges have the same overlap)
        
        This removes overlapping bases exactly once (on the parent side) for each link.
        Warnings are issued if the assumption is violated, and the max overlap is used.
        
        This modifies the node sequences in place.
        """
        # Check for violations of our assumption and collect trim values
        suffix_trims: Dict[str, int] = {}
        warning_count = 0
        
        # Process nodes with outgoing edges (parents)
        for node_id in self.all_nodes:
            if node_id not in self.nodes:
                continue
                
            outgoing = self.out_edges.get(node_id, [])
            if not outgoing:
                continue
            
            # Get all overlaps for outgoing edges
            overlaps = [overlap for _, _, _, overlap in outgoing]
            
            if overlaps:
                # Check if all overlaps are the same (our assumption)
                unique_overlaps = set(overlaps)
                if len(unique_overlaps) > 1:
                    print(f"WARNING: Node {node_id} has {len(outgoing)} outgoing edges with different overlaps: {sorted(unique_overlaps)}")
                    warning_count += 1
                    # Use max overlap to be safe
                    overlap_value = max(overlaps)
                else:
                    overlap_value = overlaps[0]
                
                if overlap_value > 0:
                    suffix_trims[node_id] = overlap_value
        
        if warning_count > 0:
            print(f"  WARNING: Found {warning_count} node(s) violating overlap assumption. Using max overlap value(s).")
        
        # Apply suffix trims (parents)
        suffix_trimmed_count = 0
        suffix_trimmed_bp = 0
        
        for node_id, trim_bp in suffix_trims.items():
            if node_id not in self.nodes:
                continue
            original_seq = self.nodes[node_id]
            if len(original_seq) > trim_bp:
                self.nodes[node_id] = original_seq[:-trim_bp]
                suffix_trimmed_count += 1
                suffix_trimmed_bp += trim_bp
        
        print(f"  Overlap resolution: trimmed {suffix_trimmed_count} nodes")
        print(f"    - Suffix trims: {suffix_trimmed_count} nodes, {suffix_trimmed_bp} bp")
        print(f"    - Total: {suffix_trimmed_bp} bp removed")
    
    def write_gfa(self, output_path: str) -> None:
        """
        Write the graph to a GFA file.
        
        Parameters
        ----------
        output_path : str
            Path to write the GFA file
        """
        with open(output_path, "w") as out:
            # Write header (optional but common)
            out.write("H\tVN:Z:1.0\n")
            
            # Write all S (segment) lines
            for node_id in sorted(self.all_nodes):
                if node_id in self.nodes:
                    sequence = self.nodes[node_id]
                    out.write(f"S\t{node_id}\t{sequence}\n")
            
            # Write all L (link) lines with overlap set to 0M
            # Use a set to avoid duplicate edges
            written_edges = set()
            for from_node in sorted(self.out_edges.keys()):
                for to_node, from_orient, to_orient, _ in self.out_edges[from_node]:
                    edge_key = (from_node, to_node, from_orient, to_orient)
                    if edge_key not in written_edges:
                        out.write(f"L\t{from_node}\t{from_orient}\t{to_node}\t{to_orient}\t0M\n")
                        written_edges.add(edge_key)

    def get_potential_start_nodes(self) -> Set[str]:
        return {n for n in self.all_nodes if len(self.in_edges[n]) == 0}

    def get_potential_end_nodes(self) -> Set[str]:
        return {n for n in self.all_nodes if len(self.out_edges[n]) == 0}

    def dfs_explore_from_start(
        self, start_node: str
    ) -> Tuple[Set[str], Dict[str, Optional[str]]]:

        visited: Set[str] = set()
        parent_map: Dict[str, Optional[str]] = {start_node: None}
        stack = [start_node]

        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)

            for to_node, _, _, _ in self.out_edges[node]:
                if to_node not in visited:
                    parent_map[to_node] = node
                    stack.append(to_node)

        return visited, parent_map

    def backtrack_path(
        self,
        start_node: str,
        end_node: str,
        parent_map: Dict[str, Optional[str]],
    ) -> Optional[List[str]]:

        if end_node not in parent_map:
            return None

        path = []
        cur = end_node
        while cur is not None:
            path.append(cur)
            cur = parent_map.get(cur)
            if len(path) > len(parent_map) + 10:
                return None

        path.reverse()
        return path if path[0] == start_node else None

    def find_all_start_to_end_paths(
        self,
        start_nodes: Optional[Set[str]] = None,
        end_nodes: Optional[Set[str]] = None,
    ) -> Dict[Tuple[str, str], List[str]]:

        if start_nodes is None:
            start_nodes = self.get_potential_start_nodes()
        if end_nodes is None:
            end_nodes = self.get_potential_end_nodes()

        paths: Dict[Tuple[str, str], List[str]] = {}

        for start in start_nodes:
            visited, parent_map = self.dfs_explore_from_start(start)
            reached_ends = visited & end_nodes

            if start in reached_ends:
                paths[(start, start)] = [start]

            for end in reached_ends:
                if end == start:
                    continue
                path = self.backtrack_path(start, end, parent_map)
                if path:
                    paths.setdefault((start, end), path)

        return paths

    def dfs_find_path_to_end_recursive(
        self,
        current_node: str,
        end_nodes: Set[str],
        visited: Set[str],
        path: List[str],
        visited_in_first_path: Set[str],
        prefer_unvisited: bool = True,
        max_depth: int = 100000000,
    ) -> Optional[List[str]]:
        """
        Recursive DFS to find a path to any end node.
        Prioritizes children not visited in first path if prefer_unvisited is True.
        
        Parameters
        ----------
        current_node : str
            Current node in the DFS
        end_nodes : Set[str]
            Set of valid end nodes
        visited : Set[str]
            Currently visited nodes in this DFS (for cycle detection)
        path : List[str]
            Current path being built
        visited_in_first_path : Set[str]
            Nodes visited in the first path
        prefer_unvisited : bool
            If True, prioritize children not in visited_in_first_path
        max_depth : int
            Maximum depth to prevent infinite recursion
            
        Returns
        -------
        Optional[List[str]]
            Complete path to end, or None if no path exists
        """
        # Prevent cycles and excessive depth
        if current_node in visited or len(path) > max_depth:
            return None
        
        visited.add(current_node)
        path.append(current_node)
        
        # Check if we reached an end
        if current_node in end_nodes:
            return path.copy()
        
        # Get children and prioritize
        children = [to_node for to_node, _, _, _ in self.out_edges[current_node]]
        if not children:
            path.pop()
            visited.remove(current_node)
            return None
        
        # Sort children: if prefer_unvisited, prioritize unvisited ones from first path
        if prefer_unvisited:
            children.sort(key=lambda x: (x in visited_in_first_path, x))
        else:
            children.sort(key=lambda x: x)
        
        # Try each child in order
        for child in children:
            result = self.dfs_find_path_to_end_recursive(
                child, end_nodes, visited, path, visited_in_first_path, prefer_unvisited, max_depth
            )
            if result:
                return result
        
        # Backtrack
        path.pop()
        visited.remove(current_node)
        return None

    def construct_second_path(
        self,
        start_node: str,
        end_nodes: Set[str],
        visited_in_first_path: Set[str],
    ) -> Optional[List[str]]:
        """
        Construct a second path that tries to avoid nodes in the first path.
        First tries with preference for unvisited nodes, then falls back to any path.
        
        Parameters
        ----------
        start_node : str
            Starting node for the second path
        end_nodes : Set[str]
            Set of valid end nodes
        visited_in_first_path : Set[str]
            Nodes visited in the first path
            
        Returns
        -------
        Optional[List[str]]
            Second path from start to end, or None if no path exists
        """
        if start_node in end_nodes:
            return [start_node]
        
        # Try with preference for unvisited nodes
        path = self.dfs_find_path_to_end_recursive(
            start_node, end_nodes, set(), [], visited_in_first_path, prefer_unvisited=True
        )
        if path:
            return path
        
        # Fallback: try without preference
        return self.dfs_find_path_to_end_recursive(
            start_node, end_nodes, set(), [], visited_in_first_path, prefer_unvisited=False
        )

    def get_path_sequence(self, path: List[str]) -> str:
        """
        Get the concatenated sequence for a path of nodes.
        
        Parameters
        ----------
        path : List[str]
            List of node IDs representing the path
            
        Returns
        -------
        str
            Concatenated sequence for the path
        """
        sequences = []
        for node_id in path:
            if node_id in self.nodes:
                sequences.append(self.nodes[node_id])
            else:
                # If node sequence is missing, use empty string
                sequences.append("")
        return "".join(sequences)


def parse_overlap(overlap_str: str) -> int:
    """
    Parse overlap string from GFA L lines.
    
    Format is typically "20M" meaning 20bp match/overlap.
    Returns the number of overlapping bases.
    """
    if not overlap_str or overlap_str == "*":
        return 0
    
    # Remove trailing 'M' if present (CIGAR-like format)
    overlap_str = overlap_str.strip()
    if overlap_str.endswith("M"):
        try:
            return int(overlap_str[:-1])
        except ValueError:
            return 0
    
    # Try parsing as plain integer
    try:
        return int(overlap_str)
    except ValueError:
        return 0


def load_gfa(gfa_path: str, resolve_overlaps: bool = True) -> GFAGraph:
    """
    Load a GFA file and optionally resolve overlaps.
    
    Parameters
    ----------
    gfa_path : str
        Path to the GFA file
    resolve_overlaps : bool
        If True, trim node sequences to resolve overlaps between linked contigs
        
    Returns
    -------
    GFAGraph
        The loaded (and optionally overlap-resolved) graph
    """
    graph = GFAGraph()

    with open(gfa_path) as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            if line.startswith("S"):
                parts = line.split("\t")
                if len(parts) >= 3:
                    node_id = parts[1]
                    # In GFA, '*' means the sequence is not present.
                    sequence_field = parts[2] if len(parts) > 2 else ""
                    sequence = "" if sequence_field == "*" else sequence_field
                    graph.add_node(node_id, sequence)
            elif line.startswith("L"):
                parts = line.split("\t")
                if len(parts) >= 6:
                    overlap = parse_overlap(parts[5])
                    graph.add_edge(parts[1], parts[2], parts[3], parts[4], overlap)

    if resolve_overlaps:
        print(f"Resolving overlaps in GFA graph...")
        graph.resolve_overlaps()

    return graph


def group_paths_into_assembled_chains(
    paths: List[Tuple[Tuple[str, str], List[str]]]
) -> Dict[int, List[Tuple[Tuple[str, str], List[str]]]]:
    """
    Group paths into assembled chains.
    Any two paths sharing at least one node belong to the same chain.
    """

    n = len(paths)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[ry] = rx

    node_to_path_indices: Dict[str, Set[int]] = defaultdict(set)
    for i, (_, path) in enumerate(paths):
        for node in path:
            node_to_path_indices[node].add(i)

    for indices in node_to_path_indices.values():
        indices = list(indices)
        for i in range(1, len(indices)):
            union(indices[0], indices[i])

    chains: Dict[int, List[Tuple[Tuple[str, str], List[str]]]] = defaultdict(list)
    for i, path in enumerate(paths):
        chains[find(i)].append(path)

    return {i: chains[k] for i, k in enumerate(sorted(chains))}


def parse_chunk_points_file(chunk_points_path: str) -> Dict[str, str]:
    chunk_to_chromosome: Dict[str, str] = {}

    with open(chunk_points_path) as f:
        f.readline()  # header
        for idx, line in enumerate(f):
            parts = line.strip().split("\t")
            if len(parts) >= 9 and "#" in parts[8]:
                chrom = parts[8].split("#")[2]
                chunk_to_chromosome[str(idx)] = chrom

    return chunk_to_chromosome


def extract_chunk_id_from_node(node_id: str) -> Optional[str]:
    return node_id.split("#")[0] if "#" in node_id else None


def group_paths_by_chromosome_and_assembled_chains(
    paths: List[Tuple[Tuple[str, str], List[str]]],
    chunk_to_chromosome: Dict[str, str],
) -> Dict[str, Dict[int, List[Tuple[Tuple[str, str], List[str]]]]]:

    paths_by_chromosome: Dict[str, List[Tuple[Tuple[str, str], List[str]]]] = defaultdict(list)

    for (start, end), path in paths:
        start_chunk = extract_chunk_id_from_node(start)
        end_chunk = extract_chunk_id_from_node(end)

        if start_chunk not in chunk_to_chromosome or end_chunk not in chunk_to_chromosome:
            raise ValueError(f"Missing chunk mapping for path {path}")

        chrom = chunk_to_chromosome[start_chunk]
        if chrom != chunk_to_chromosome[end_chunk]:
            raise ValueError(f"Cross-chromosome path: {path}")

        paths_by_chromosome[chrom].append(((start, end), path))

    # return format: {chrom: {assembled_chain_id: [(start, end), path], ...}, ...}
    return {
        chrom: group_paths_into_assembled_chains(plist)
        for chrom, plist in paths_by_chromosome.items()
    }


def construct_sibling_dicts(graph: GFAGraph, contigs: List[str], is_start: bool = True) -> Dict[str, List[str]]:
    """
    Try to classify siblings
    Assembled segment names are of the form a-b-c-d-Pn, where:
        a-b identifies the bubble chain.
        c is the position of the bubble in the bubble chain.
        d identifies the haplotype in the bubble.
        n is the ploidy of the bubble.
    """

    sibling_dict = {k: [] for k in contigs}

    # Compare all pairs of contigs
    for i in range(len(contigs)):
        a_contig = contigs[i]
        a_parts = a_contig.split('#')[1].split('-')
        
        for j in range(i + 1, len(contigs)):
            other_contig = contigs[j]
            o_parts = other_contig.split('#')[1].split('-')
            
            if a_parts[4] not in ["P0", "P1"]:
                # if everything is same, except part[3], then they are siblings
                if a_parts[0:3] == o_parts[0:3] and a_parts[4] == o_parts[4] and a_parts[3] != o_parts[3]:
                    sibling_dict[a_contig].append(other_contig)
                    sibling_dict[other_contig].append(a_contig)
            
            elif a_parts[4] == "P0":
                # if everything is same, except part[1], then they are siblings. Note, they have to be incremental
                if a_parts[0] == o_parts[0] and a_parts[2:5] == o_parts[2:5] and (
                    int(a_parts[1]) + 1 == int(o_parts[1]) or
                    int(a_parts[1]) - 1 == int(o_parts[1])
                ):
                    # Check if they are connected to atleast one same contig
                    if is_start:
                        # check if they share atleast one child contig
                        a_contig_children = [child for (child, _, _, _) in graph.out_edges[a_contig]]
                        o_contig_children = [child for (child, _, _, _) in graph.out_edges[other_contig]]
                        if len(set(a_contig_children) & set(o_contig_children)) > 0:
                            sibling_dict[a_contig].append(other_contig)
                            sibling_dict[other_contig].append(a_contig)
                    else:
                        # check if they share atleast one parent contig
                        a_contig_parents = [parent for (parent, _, _, _) in graph.in_edges[a_contig]]
                        o_contig_parents = [parent for (parent, _, _, _) in graph.in_edges[other_contig]]
                        if len(set(a_contig_parents) & set(o_contig_parents)) > 0:
                            sibling_dict[a_contig].append(other_contig)
                            sibling_dict[other_contig].append(a_contig)

            # elif a_parts[4] == "P0":
            #     # if everything is same, except part[1], then they are siblings. Note, they have to be incremental
            #     if a_parts[0] == o_parts[0] and a_parts[2:5] == o_parts[2:5] and (
            #         int(a_parts[1]) + 1 == int(o_parts[1]) or
            #         int(a_parts[1]) - 1 == int(o_parts[1])
            #     ):
            #         sibling_dict[a_contig].append(other_contig)
            #         sibling_dict[other_contig].append(a_contig)
    
    return sibling_dict


def select_longest_paths_per_assembled_chain_candidate(graph: GFAGraph, chromosome_assembled_chain_candidates: Dict[str, Dict[int, List[Tuple[Tuple[str, str], List[str]]]]]):
    """
    Select the longest path per assembled chain candidate.
    If there are multiple paths with the same length, select all of them.
    """
    longest_paths_per_assembled_chain_candidate = defaultdict(dict)

    for chrom, chains in chromosome_assembled_chain_candidates.items():
        for cid, plist in chains.items():
            # sort the paths by length, and select the longest path
            sorted_paths = sorted(plist, key=lambda x: len(x[1]), reverse=True)
            # select longest path's length.
            longest_path_length = len(sorted_paths[0][1])
            # check if there are multiple paths (could be at max 4 in theory) with the same length, and if so, select any one of them, 
            # but report a WARNING
            paths_with_longest_path_length = [p for p in sorted_paths if len(p[1]) == longest_path_length]
            # Track already added paths to avoid duplicates
            added_paths = {(p[0][0], p[0][1]) for p in paths_with_longest_path_length}

            # Check if we missed a sibling
            # Keep only unique starts and ends
            starts = list(set([p[0][0] for p in sorted_paths]))
            ends = list(set([p[0][1] for p in sorted_paths]))
            # Construct sibling dicts for starts and ends
            starts_sibling_dict = construct_sibling_dicts(graph, starts, is_start=True)
            ends_sibling_dict = construct_sibling_dicts(graph, ends, is_start=False)

            # Using sibling dicts, if only one sibling is found in the longest path(s), 
            # then select the path with the other sibling as well.
            starts_in_longest_paths = [p[0][0] for p in paths_with_longest_path_length]
            ends_in_longest_paths = [p[0][1] for p in paths_with_longest_path_length]
            starts_in_longest_paths_sibling_dict = construct_sibling_dicts(graph, starts_in_longest_paths, is_start=True)
            ends_in_longest_paths_sibling_dict = construct_sibling_dicts(graph, ends_in_longest_paths, is_start=False)
            
            for start in starts_in_longest_paths:
                if len(starts_in_longest_paths_sibling_dict[start]) < len(starts_sibling_dict[start]):
                    print(f"FOUND MISSING START SIBLING FOR {start}")
                    # it means we missed a sibling
                    # so, we need to add the missing sibling's longest path to the longest path(s)
                    missing_siblings = [s for s in starts_sibling_dict[start] if s not in starts_in_longest_paths]
                    print(f"  Missing siblings for start={start} (not in starts_in_longest_paths): {missing_siblings}")
                    for sibling in missing_siblings:
                        sibling_paths = [p for p in sorted_paths if p[0][0] == sibling]
                        if sibling_paths:
                            sibling_path = sibling_paths[0]
                            path_key = (sibling_path[0][0], sibling_path[0][1])
                            if path_key not in added_paths:
                                print(f"  ADDING sibling path: {path_key}")
                                paths_with_longest_path_length.append(sibling_path)
                                added_paths.add(path_key)
                            else:
                                print(f"  SKIPPING (already exists): {path_key}")
                        else:
                            print(f"  No paths found for start sibling: {sibling}")
            
            for end in ends_in_longest_paths:
                if len(ends_in_longest_paths_sibling_dict[end]) < len(ends_sibling_dict[end]):
                    print(f"FOUND MISSING END SIBLING FOR {end}")
                    missing_siblings = [e for e in ends_sibling_dict[end] if e not in ends_in_longest_paths]
                    print(f"  Missing siblings for end={end} (not in ends_in_longest_paths): {missing_siblings}")
                    for sibling in missing_siblings:
                        sibling_paths = [p for p in sorted_paths if p[0][1] == sibling]
                        if sibling_paths:
                            sibling_path = sibling_paths[0]
                            path_key = (sibling_path[0][0], sibling_path[0][1])
                            if path_key not in added_paths:
                                print(f"  ADDING sibling path: {path_key}")
                                paths_with_longest_path_length.append(sibling_path)
                                added_paths.add(path_key)
                            else:
                                print(f"  SKIPPING (already exists): {path_key}")
                        else:
                            print(f"  No paths found for end sibling: {sibling}")

            longest_paths_per_assembled_chain_candidate[chrom][cid] = paths_with_longest_path_length

    return longest_paths_per_assembled_chain_candidate


def construct_pseudohaplotypes(
    graph: GFAGraph,
    chain_bounds: Dict[str, Dict[int, Tuple[List[str], List[str]]]],
) -> Dict[str, Dict[int, Tuple[List[str], List[str]]]]:
    """
    Construct pseudohaplotypes for each assembled chain.
    For each chain, construct 2 haplotypes by finding complementary paths.
    
    Parameters
    ----------
    graph : GFAGraph
        The graph structure
    chain_bounds : Dict[str, Dict[int, Tuple[List[str], List[str]]]]
        Dictionary mapping chromosome -> chain_id -> (starts, ends)
        
    Returns
    -------
    Dict[str, Dict[int, Tuple[List[str], List[str]]]]
        Dictionary mapping chromosome -> chain_id -> (haplotype1, haplotype2)
    """
    pseudohaplotypes: Dict[str, Dict[int, Tuple[List[str], List[str]]]] = defaultdict(dict)
    
    for chrom, chains in chain_bounds.items():
        for cid, (starts, ends) in chains.items():
            starts_set = set(starts)
            ends_set = set(ends)
            
            if not starts_set or not ends_set:
                print(f"Warning: Chain {chrom}:{cid} has no starts or ends, skipping")
                continue
            
            # Step 1: Find first path from any start to any end
            first_path = None
            first_start = None
            first_end = None
            
            for start in starts_set:
                # Try simple DFS to any end
                for end in ends_set:
                    visited, parent_map = graph.dfs_explore_from_start(start)
                    if end in visited:
                        first_path = graph.backtrack_path(start, end, parent_map)
                        if first_path:
                            first_start = start
                            first_end = end
                            break
                if first_path:
                    break
            
            if not first_path:
                print(f"Warning: Could not find first path for chain {chrom}:{cid}, skipping")
                continue
            
            visited_in_first_path = set(first_path)
            
            # Step 2: Find second path
            # Pick second start (or same if only one)
            second_start = None
            if len(starts_set) > 1:
                # Pick a different start
                for start in starts_set:
                    if start != first_start:
                        second_start = start
                        break
            else:
                # Use same start if only one available
                second_start = first_start
            
            # Construct second path
            second_path = graph.construct_second_path(
                second_start, ends_set, visited_in_first_path
            )
            
            if not second_path:
                print(f"Warning: Could not find second path for chain {chrom}:{cid}, using first path for both")
                second_path = first_path.copy()
            
            # Store the two haplotypes
            pseudohaplotypes[chrom][cid] = (first_path, second_path)
    
    return pseudohaplotypes


def compute_ploidy_assignments(
    pseudohaplotypes: Dict[str, Dict[int, Tuple[List[str], List[str]]]],
    graph: GFAGraph,
) -> List[Tuple[str, str, str, str, str]]:
    """
    Compute ploidy assignments for contigs based on pseudohaplotype construction.
    A contig is haploid if it appears in both hap1 and hap2 for the same chain
    (same contig chosen in both paths, indicating only one version exists).
    A contig is diploid if it appears in only one haplotype
    (different contigs chosen in different haplotypes, indicating two versions exist).
    
    Parameters
    ----------
    pseudohaplotypes : Dict[str, Dict[int, Tuple[List[str], List[str]]]]
        Dictionary mapping chromosome -> chain_id -> (haplotype1, haplotype2)
    graph : GFAGraph
        The graph structure (for checking if nodes have alternative children)
        
    Returns
    -------
    List[Tuple[str, str, str, str, str, int]]
        List of (chromosome, chain_id, contig_id, ploidy, context, length) tuples
        where ploidy is "haploid" or "diploid", context describes why it's classified that way,
        and length is the sequence length of the contig
    """
    ploidy_assignments = []
    
    for chrom, chains in sorted(pseudohaplotypes.items()):
        for cid, (hap1, hap2) in sorted(chains.items()):
            # Get sets of nodes in each haplotype
            hap1_nodes = set(hap1)
            hap2_nodes = set(hap2)
            
            # All nodes that appear in either haplotype
            all_nodes = hap1_nodes | hap2_nodes
            
            # Nodes that appear in both haplotypes (haploid - same contig in both)
            haploid_nodes = hap1_nodes & hap2_nodes
            
            # Nodes that appear in only one haplotype (diploid - different contigs)
            diploid_nodes = all_nodes - haploid_nodes
            
            # Classify each node
            for node in sorted(all_nodes):
                # Get contig length
                contig_length = len(graph.nodes.get(node, ""))
                
                if node in haploid_nodes:
                    # Haploid: appears in both haplotypes (same contig chosen in both paths)
                    context_parts = []
                    
                    # Get parents of this node that appear in the paths
                    parents = [p for p, _, _, _ in graph.in_edges[node] if p in all_nodes]
                    
                    for parent in parents:
                        # Get all children of this parent
                        children = [c for c, _, _, _ in graph.out_edges[parent]]
                        if len(children) > 1:
                            # Parent had multiple children - check which ones were used
                            children_in_hap1 = [c for c in children if c in hap1_nodes]
                            children_in_hap2 = [c for c in children if c in hap2_nodes]
                            
                            # If parent had alternatives but same child was chosen in both paths
                            if node in children_in_hap1 and node in children_in_hap2:
                                alternatives = [c for c in children if c != node]
                                context_parts.append(
                                    f"parent_{parent}_had_{len(alternatives)}_alternatives_but_same_chosen"
                                )
                    
                    if not context_parts:
                        # Still haploid, but no parent context available
                        context = "appears_in_both_haplotypes"
                    else:
                        context = "; ".join(context_parts)
                    
                    ploidy_assignments.append((chrom, str(cid), node, "haploid", context, contig_length))
                else:
                    # Diploid: appears in only one haplotype (different contigs in different haplotypes)
                    which_hap = "hap1" if node in hap1_nodes else "hap2"
                    ploidy_assignments.append((chrom, str(cid), node, "diploid", f"appears_only_in_{which_hap}", contig_length))
    
    return ploidy_assignments


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("gfa")
    parser.add_argument("-c", "--chunk-points", required=True)
    parser.add_argument("-o", "--output-prefix", required=True)
    args = parser.parse_args()

    print("\nLoading GFA file and resolving overlaps...")
    t0 = time.time()
    graph = load_gfa(args.gfa)
    elapsed = time.time() - t0
    print(f"GFA loading and overlap resolution finished in {elapsed:.2f} seconds")
    
    # Write overlap-resolved GFA file
    print(f"\nWriting overlap-resolved GFA file...")
    t0 = time.time()
    gfa_dir = os.path.dirname(args.gfa)
    gfa_basename = os.path.basename(args.gfa)
    gfa_name_without_ext = os.path.splitext(gfa_basename)[0]
    resolved_gfa_path = os.path.join(gfa_dir, f"{gfa_name_without_ext}_ovlp_resolved.gfa")
    graph.write_gfa(resolved_gfa_path)
    elapsed = time.time() - t0
    print(f"Writing overlap-resolved GFA finished in {elapsed:.2f} seconds")
    print(f"  Output: {resolved_gfa_path}")
    
    print("\nFinding all start-to-end paths...")
    t0 = time.time()
    start_nodes = graph.get_potential_start_nodes()
    end_nodes = graph.get_potential_end_nodes()
    paths_dict = graph.find_all_start_to_end_paths(start_nodes, end_nodes)
    paths = list(paths_dict.items())
    elapsed = time.time() - t0
    print(f"Path finding finished in {elapsed:.2f} seconds")
    print(f"  Found {len(paths)} paths from {len(start_nodes)} start nodes to {len(end_nodes)} end nodes")

    print("\nParsing chunk points file...")
    t0 = time.time()
    chunk_to_chromosome = parse_chunk_points_file(args.chunk_points)
    elapsed = time.time() - t0
    print(f"Chunk points parsing finished in {elapsed:.2f} seconds")
    print(f"  Mapped {len(chunk_to_chromosome)} chunks to chromosomes")

    print("\nGrouping paths into assembled chain candidates...")
    t0 = time.time()
    chromosome_assembled_chain_candidates = group_paths_by_chromosome_and_assembled_chains(
        paths, chunk_to_chromosome
    )
    elapsed = time.time() - t0
    print(f"Assembled chain candidate grouping finished in {elapsed:.2f} seconds")

    total_assembled_chain_candidates = sum(len(v) for v in chromosome_assembled_chain_candidates.values())
    print(f"  Found {total_assembled_chain_candidates} assembled chain candidates")

    print("\nWriting assembled chain candidates...")
    t0 = time.time()
    with open(f"{args.output_prefix}_assembled_chain_candidates.tsv", "w") as out:
        out.write("chromosome\tassembled_chain_id\tstart\tend\tpath\n")
        for chrom, chains in chromosome_assembled_chain_candidates.items():
            for cid, plist in chains.items():
                for (s, e), p in plist:
                    out.write(f"{chrom}\t{cid}\t{s}\t{e}\t{','.join(p)}\n")
    elapsed = time.time() - t0
    print(f"Writing assembled chain candidates finished in {elapsed:.2f} seconds")

    print("\nSelecting longest paths per assembled chain candidate...")
    t0 = time.time()
    assembled_chains = select_longest_paths_per_assembled_chain_candidate(graph, chromosome_assembled_chain_candidates)
    elapsed = time.time() - t0
    print(f"Longest path selection finished in {elapsed:.2f} seconds")
    
    total_assembled_chains = sum(len(v) for v in assembled_chains.values())
    print(f"  Selected {total_assembled_chains} assembled chains")
    
    print("\nWriting assembled chains...")
    t0 = time.time()
    with open(f"{args.output_prefix}_assembled_chains.tsv", "w") as out:
        out.write("chromosome\tassembled_chain_id\tstart\tend\tpath\n")
        for chrom, chains in assembled_chains.items():
            for cid, plist in chains.items():
                for (start, end), path in plist:
                    out.write(f"{chrom}\t{cid}\t{start}\t{end}\t{','.join(path)}\n")
    elapsed = time.time() - t0
    print(f"Writing assembled chains finished in {elapsed:.2f} seconds")

    print("\nExtracting assembled chain bounds...")
    t0 = time.time()
    chain_bounds = defaultdict(dict)
    for chrom, chains in assembled_chains.items():
        for cid, plist in chains.items():
            starts = []
            ends = []
            for (start, end), path in plist:
                starts.append(start)
                ends.append(end)
            chain_bounds[chrom][cid] = (list(set(starts)), list(set(ends)))
    elapsed = time.time() - t0
    print(f"Chain bounds extraction finished in {elapsed:.2f} seconds")
    
    print("\nWriting assembled chain bounds...")
    t0 = time.time()
    with open(f"{args.output_prefix}_assembled_chain_bounds.tsv", "w") as out:
        out.write("chromosome\tassembled_chain_id\tstart\tend\n")
        for chrom, chains in chain_bounds.items():
            for cid, (starts, ends) in chains.items():
                out.write(f"{chrom}\t{cid}\t{','.join(starts)}\t{','.join(ends)}\n")
    elapsed = time.time() - t0
    print(f"Writing chain bounds finished in {elapsed:.2f} seconds")


    ## Construct pseudo-haplotypes
    print(f"\nConstructing pseudohaplotypes for {len(chain_bounds)} chromosomes...")
    t0 = time.time()
    pseudohaplotypes = construct_pseudohaplotypes(graph, chain_bounds)
    elapsed = time.time() - t0
    
    total_pseudohaplotypes = sum(len(chains) for chains in pseudohaplotypes.values())
    print(f"Pseudohaplotype construction finished in {elapsed:.2f} seconds")
    print(f"  Successfully constructed {total_pseudohaplotypes} pseudohaplotype pairs")
    
    # Write pseudohaplotypes to output file
    print("\nWriting pseudohaplotypes TSV file...")
    t0 = time.time()
    pseudohaplotypes_output_file = f"{args.output_prefix}_pseudohaplotypes.tsv"
    with open(pseudohaplotypes_output_file, "w") as out:
        out.write("chromosome\tassembled_chain_id\thaplotype_id\tstart\tend\tpath\n")
        for chrom, chains in sorted(pseudohaplotypes.items()):
            for cid, (hap1, hap2) in sorted(chains.items()):
                out.write(f"{chrom}\t{cid}\t1\t{hap1[0]}\t{hap1[-1]}\t{','.join(hap1)}\n")
                out.write(f"{chrom}\t{cid}\t2\t{hap2[0]}\t{hap2[-1]}\t{','.join(hap2)}\n")
    elapsed = time.time() - t0
    print(f"Writing pseudohaplotypes TSV finished in {elapsed:.2f} seconds")
    print(f"  Output: {pseudohaplotypes_output_file}")
    
    # Write FASTA files for hap1 and hap2
    print("\nWriting FASTA files for hap1 and hap2...")
    t0 = time.time()
    hap1_fasta_file = f"{args.output_prefix}_hap1.fasta"
    hap2_fasta_file = f"{args.output_prefix}_hap2.fasta"
    
    with open(hap1_fasta_file, "w") as hap1_out, open(hap2_fasta_file, "w") as hap2_out:
        for chrom, chains in sorted(pseudohaplotypes.items()):
            for cid, (hap1, hap2) in sorted(chains.items()):
                # Get sequences for hap1
                hap1_seq = graph.get_path_sequence(hap1)
                hap1_header = f">{chrom}_chain_{cid}_hap1_start_{hap1[0]}_end_{hap1[-1]}"
                hap1_out.write(f"{hap1_header}\n")
                # Write sequence in 80 character lines (standard FASTA format)
                for i in range(0, len(hap1_seq), 80):
                    hap1_out.write(f"{hap1_seq[i:i+80]}\n")
                
                # Get sequences for hap2
                hap2_seq = graph.get_path_sequence(hap2)
                hap2_header = f">{chrom}_chain_{cid}_hap2_start_{hap2[0]}_end_{hap2[-1]}"
                hap2_out.write(f"{hap2_header}\n")
                # Write sequence in 80 character lines (standard FASTA format)
                for i in range(0, len(hap2_seq), 80):
                    hap2_out.write(f"{hap2_seq[i:i+80]}\n")
    elapsed = time.time() - t0
    print(f"Writing FASTA files finished in {elapsed:.2f} seconds")
    print(f"  Output: {hap1_fasta_file}, {hap2_fasta_file}")
    
    # Compute and write ploidy assignments
    print(f"\nComputing ploidy assignments...")
    t0 = time.time()
    ploidy_assignments = compute_ploidy_assignments(pseudohaplotypes, graph)
    elapsed = time.time() - t0
    
    diploid_count = sum(1 for _, _, _, ploidy, _, _ in ploidy_assignments if ploidy == "diploid")
    haploid_count = sum(1 for _, _, _, ploidy, _, _ in ploidy_assignments if ploidy == "haploid")
    print(f"Ploidy assignment computation finished in {elapsed:.2f} seconds")
    print(f"  Found {diploid_count} diploid contigs and {haploid_count} haploid contigs")
    
    print("\nWriting ploidy assignments TSV file...")
    t0 = time.time()
    ploidy_output_file = f"{args.output_prefix}_ploidy_assignments.tsv"
    with open(ploidy_output_file, "w") as out:
        out.write("chromosome\tassembled_chain_id\tcontig_id\tploidy\tcontext\tlength\n")
        for chrom, chain_id, contig_id, ploidy, context, length in ploidy_assignments:
            out.write(f"{chrom}\t{chain_id}\t{contig_id}\t{ploidy}\t{context}\t{length}\n")
    elapsed = time.time() - t0
    print(f"Writing ploidy assignments TSV finished in {elapsed:.2f} seconds")
    print(f"  Output: {ploidy_output_file}")
    
    print("\nWriting Bandage CSV files...")
    t0 = time.time()
    # Create ploidy mapping: contig_id -> ploidy
    contig_to_ploidy = {}
    for chrom, chain_id, contig_id, ploidy, context, length in ploidy_assignments:
        # If a contig appears in multiple chains, use the first ploidy assignment
        # (in practice, a contig should have the same ploidy across chains)
        if contig_id not in contig_to_ploidy:
            contig_to_ploidy[contig_id] = ploidy
    
    # Collect all contigs from hap1 and hap2 paths
    hap1_contigs = set()
    hap2_contigs = set()
    for chrom, chains in pseudohaplotypes.items():
        for cid, (hap1, hap2) in chains.items():
            hap1_contigs.update(hap1)
            hap2_contigs.update(hap2)
    
    # Write Bandage CSV files
    hap1_bandage_file = f"{args.output_prefix}_hap1_contigs.bandage.csv"
    hap2_bandage_file = f"{args.output_prefix}_hap2_contigs.bandage.csv"
    
    with open(hap1_bandage_file, "w") as out:
        out.write("Contig,Ploidy,Color\n")
        for contig in sorted(hap1_contigs):
            ploidy = contig_to_ploidy.get(contig, "unknown")
            out.write(f"{contig},{ploidy},#bf3030\n")
    
    with open(hap2_bandage_file, "w") as out:
        out.write("Contig,Ploidy,Color\n")
        for contig in sorted(hap2_contigs):
            ploidy = contig_to_ploidy.get(contig, "unknown")
            out.write(f"{contig},{ploidy},#30bf30\n")
    elapsed = time.time() - t0
    print(f"Writing Bandage CSV files finished in {elapsed:.2f} seconds")
    print(f"  Output: {hap1_bandage_file}, {hap2_bandage_file}")


if __name__ == "__main__":
    main()
