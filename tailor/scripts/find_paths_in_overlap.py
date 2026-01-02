#!/usr/bin/env python3
"""
find_paths_in_overlap.py

Steps:
- From PAF, finds the longest alignment for each query-path:ref-path combination
- Finds the alignment with maximum matches across all selected combinations -> this is the Hap1 path
- Now, we found one path (which will represent one haplotype), and need another path, 
    which will represent Hap-2. So, for the other haplotype, we should select the path(s) 
    with corresponding sibling contigs from the query path. To find best such path, 
    choose the one that has the maximum matches.

Output: Print
"""

import argparse
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass


@dataclass
class Alignment:
    """A single alignment from PAF."""
    q_name: str       # Query path name
    q_len: int
    q_start: int
    q_end: int
    strand: str
    t_name: str       # Target/ref path name
    t_len: int
    t_start: int
    t_end: int
    num_matches: int
    aln_block_len: int
    mapq: int
    cs_line: List[Tuple[str, int]]
    line: str

@dataclass
class Contig:
    """A single contig from GFA."""
    name: str
    chunk_id: str
    length: int
    sequence: str
    links: List[Tuple[str, str, int]]   # (to_contig, orientation, overlap)
    
    def __repr__(self):
        return f"Contig(name={self.name}, chunk_id={self.chunk_id}, length={self.length}), links={self.links}"


def read_gfa_file_for_candidate_contigs(gfa_path: str, candidate_contig_names: Set[str]) -> Dict[str, Contig]:
    """Read GFA file and return dictionary of contigs."""
    contigs = {}

    with open(gfa_path, 'r') as f:
        for line in f:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                # print(f"  At contig: {parts[1]}")
                if parts[1] in candidate_contig_names:  # Only extract info for candidate contigs
                    contigs[parts[1]] = Contig(name=parts[1], chunk_id="query", length=int(parts[3].split("LN:i:")[-1]), sequence=parts[2], links=[])
            
            if line.startswith('L'):    # Assuming correct GFA strutcure that L lines follow S lines
                link_parts = line.strip().split('\t')
                if link_parts[1] in candidate_contig_names:
                    contigs[link_parts[1]].links.append((link_parts[3], link_parts[4], int(link_parts[5][:-1])))
    
    return contigs


def get_cs_line(cigar: str) -> List[Tuple[str, int]]:
    """Get list of (operation, length) tuples from CIGAR string."""
    cs_string = cigar.split('cg:Z:')[1]
    cs_line = []
    count_str = ""
    for i in cs_string:
        if i.isdigit():
            count_str += i
        else:
            if count_str:  # Only append if we have a count
                cs_line.append((i, int(count_str)))
                count_str = ""
    return cs_line

def parse_paf(paf_path: str) -> List[Alignment]:
    """Parse PAF file and return list of Alignment objects."""
    alignments = []
    with open(paf_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')       
            try:
                aln = Alignment(
                    q_name=parts[0],
                    q_len=int(parts[1]),
                    q_start=int(parts[2]),
                    q_end=int(parts[3]),
                    strand=parts[4],
                    t_name=parts[5],
                    t_len=int(parts[6]),
                    t_start=int(parts[7]),
                    t_end=int(parts[8]),
                    num_matches=int(parts[9]),
                    aln_block_len=int(parts[10]),
                    mapq=int(parts[11]) if len(parts) > 11 else 0,
                    cs_line=get_cs_line(parts[-1]),
                    line=line,
                )
                alignments.append(aln)
            except (ValueError, IndexError) as e:
                print(f"Warning: Skipping malformed PAF line: {e}")
                continue
    
    return alignments

def check_if_contigs_are_siblings(contig1_name: str, contig2_name: str) -> bool:
    contig1_list = contig1_name.split('-')
    contig2_list = contig2_name.split('-')
    # our mental model for checking siblings is very heuristic and error-prone.
    # for now, we believe sibling contigs should have values at indices 1 and 2 be the same (0-indexed), and the value at index 3 be different (which represents the phase according to shasta).
    if len(contig1_list) != len(contig2_list) or (len(contig1_list) != 5):
        return False
    if (contig1_list[1] != contig2_list[1]) or (contig1_list[2] != contig2_list[2]):
        return False
    if (contig1_list[3] != contig2_list[3]):
        return True
    return False


def create_grouped_sibling_contigs_dict(connectivity_dict: Dict[int, List[str]]) -> Dict[int, List[List[str]]]:
    sibling_grouped_connectivity_dict = {}
    for idx, contig_list in connectivity_dict.items():
        connected_component_idx = 0
        sibling_contig_map = {} # key: contig name, value: index of connected component to which it belongs
        contig_name_and_idx_list = []
        for contig_name in contig_list:
            for contig_name2 in contig_list:
                if contig_name == contig_name2:
                    break
                if check_if_contigs_are_siblings(contig_name, contig_name2):
                    sibling_contig_map[contig_name] = sibling_contig_map[contig_name2]
                    contig_name_and_idx_list.append((contig_name, sibling_contig_map[contig_name]))
                    break
            if contig_name not in sibling_contig_map:
                sibling_contig_map[contig_name] = connected_component_idx
                contig_name_and_idx_list.append((contig_name, sibling_contig_map[contig_name]))
                connected_component_idx += 1
        
        # now we convert the sibling_contig_map to a list of lists of sibling contigs
        sorted_contig_name_and_idx_list = sorted(contig_name_and_idx_list, key=lambda x: x[1])
        sibling_grouped_connectivity_dict[idx] = []
        curr_component = []
        curr_component_idx = 0
        for contig_name, contig_idx in sorted_contig_name_and_idx_list:
            if contig_idx != curr_component_idx:
                sibling_grouped_connectivity_dict[idx].append(curr_component)
                curr_component_idx = contig_idx
                curr_component = []
            curr_component.append(contig_name)
        sibling_grouped_connectivity_dict[idx].append(curr_component)
    return sibling_grouped_connectivity_dict

def generate_contig_relationships_dict(alignments: List[Alignment]) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Dict[int, List[List[str]]], Dict[int, List[List[str]]]]:
    """
    Generate a dictionary of contig relationships from the alignments.
    For example: 
        {0: [0-8-2-0-P1], 1: [0-8-3-0-P3, 0-8-3-1-P3, 0-8-3-2-P3], 2: [0-8-3-0-P3, 0-8-3-2-P3], ...}
        This means that contig 0-8-2-0-P1 is a haploid contig, and
        0-8-3-0-P3, 0-8-3-1-P3, 0-8-3-2-P3 are all sibling contigs
    """
    contig_relationships = {}
    # list all unique query and reference path names
    query_paths = list(set(aln.q_name for aln in alignments))
    ref_paths = list(set(aln.t_name for aln in alignments))

    # FOR QUERY PATHS:
    first_path = query_paths[0]
    prefix, suffix = first_path.split(':')[0], first_path.split(':')[1]
    query_connectivity_dict = {idx: [contig_name] for idx, contig_name in enumerate(suffix.split('_'))}

    for path_name in query_paths[1:]:
        prefix, suffix = path_name.split(':')[0], path_name.split(':')[1]
        contig_names = suffix.split('_')
        for idx, contig_name in enumerate(contig_names):
            if idx not in query_connectivity_dict:
                query_connectivity_dict[idx] = []
            if contig_name not in query_connectivity_dict[idx]:   # it means this is a sibling contig
                query_connectivity_dict[idx].append(contig_name)
    # Now, at each index, we group the true sibling contigs together. Thus our resulting sibling_grouped_query_connectivity dict will have the pos idx as keys, and a list of lists of sibling contigs as values.
    # Note that here we are having to check for true siblingness because query_connectivity_dict is created from the paths in the alignments, which might have differing num contigs for each path, and so the indices we assigned to them might be wrong.
    # We encountered this issue in chunk 74_75. The problem was that a haploid contig 75#0-10-2-0-P1 was wrongly being assigned 2 siblings ('75#0-11-0-0-P2', '75#0-11-0-1-P2'). This was causing that haploid contig to be skipped when creating hap2 possible qpaths (because it was being considered as a non-haploid contig, for which we create all paths maximally different from the selected hap1 path). 
    sibling_grouped_query_connectivity_dict = create_grouped_sibling_contigs_dict(query_connectivity_dict)
    print(f"  Sibling grouped query connectivity dict: \n\t{sibling_grouped_query_connectivity_dict}")
    # FOR REFERENCE PATHS:
    first_path = ref_paths[0]
    prefix, suffix = first_path.split(':')[0], first_path.split(':')[1]
    ref_connectivity_dict = {idx: [contig_name] for idx, contig_name in enumerate(suffix.split('_'))}

    for path_name in ref_paths[1:]:
        prefix, suffix = path_name.split(':')[0], path_name.split(':')[1]
        contig_names = suffix.split('_')
        for idx, contig_name in enumerate(contig_names):
            if idx not in ref_connectivity_dict:
                ref_connectivity_dict[idx] = []
            if contig_name not in ref_connectivity_dict[idx]:   # it means this is a sibling contig
                ref_connectivity_dict[idx].append(contig_name)
    sibling_grouped_ref_connectivity_dict = create_grouped_sibling_contigs_dict(ref_connectivity_dict)
    print(f"  Sibling grouped ref connectivity dict: \n\t{sibling_grouped_ref_connectivity_dict}")

    return query_connectivity_dict, ref_connectivity_dict, sibling_grouped_query_connectivity_dict, sibling_grouped_ref_connectivity_dict


def construct_hap2_path_recursively(sibling_grouped_connectivity_dict, hap1_name_contig_names, hap2_name_possible_suffixes, current_hap2_path, idx):
    # earlier we were comparing idx with len(hap1_name_contig_names) to figure out if end of hap2 path was reached. 
    # The assumption there was that hap1 and hap2 paths will be equally long. But that's wrong (for reference, check 55_56) 
    # So, now, we append ALL possible suffixes to the hap2_name_possible_suffixes list.
    if len(hap2_name_possible_suffixes) > 10000:
        raise ValueError(f"Hap2 path possible suffixes list is too long (while backtracking). This is a potential centromere region. Investigate!!!")
    if idx > max(sibling_grouped_connectivity_dict.keys()):
        hap2_name_possible_suffixes.append('_'.join(current_hap2_path))
        return
    if len(current_hap2_path) > 0:
        hap2_name_possible_suffixes.append('_'.join(current_hap2_path))
    for sibling_group in sibling_grouped_connectivity_dict[idx]:
        if len(sibling_group) == 1:
            current_hap2_path.append(sibling_group[0])
            construct_hap2_path_recursively(sibling_grouped_connectivity_dict, hap1_name_contig_names, hap2_name_possible_suffixes, current_hap2_path, idx + 1)
            current_hap2_path.pop()
        else:
            for contig_name in sibling_group:
                if (idx >= len(hap1_name_contig_names)) or ((idx < len(hap1_name_contig_names)) and (contig_name != hap1_name_contig_names[idx])):
                    current_hap2_path.append(contig_name)
                    construct_hap2_path_recursively(sibling_grouped_connectivity_dict, hap1_name_contig_names, hap2_name_possible_suffixes, current_hap2_path, idx + 1)
                    current_hap2_path.pop()
    return


def construct_hap2_paths_helper(sibling_grouped_connectivity_dict, hap1_name_contig_names, hap2_name_possible_suffixes):
    """
    Helper function to construct the hap2 qname and tname paths using the connectivity dictionary.
    """
    current_hap2_path = []
    construct_hap2_path_recursively(sibling_grouped_connectivity_dict, hap1_name_contig_names, hap2_name_possible_suffixes, current_hap2_path, 0)
    return


def construct_hap2_paths(sibling_grouped_query_connectivity_dict, sibling_grouped_ref_connectivity_dict, hap1_aln_record) -> Tuple[List[str], List[str]]:
    """
    Construct the hap2 qname and tname paths using the connectivity dictionary.
    Find the best alignment amongst selected path pairs (qname:tname) --> Hap2 path
    """
    hap1_qname_contig_names = hap1_aln_record.q_name.split(':')[1].split('_')
    hap1_tname_contig_names = hap1_aln_record.t_name.split(':')[1].split('_')

    # FOR QUERY PATH
    hap2_qname_possible_suffixes = []
    construct_hap2_paths_helper(sibling_grouped_query_connectivity_dict, hap1_qname_contig_names, hap2_qname_possible_suffixes)
    print(f"  Hap2 qname possible suffixes: \n\t{hap2_qname_possible_suffixes}")

    # FOR REFERENCE PATH
    hap2_tname_possible_suffixes = []
    construct_hap2_paths_helper(sibling_grouped_ref_connectivity_dict, hap1_tname_contig_names, hap2_tname_possible_suffixes)
    print(f"  Hap2 tname possible suffixes: \n\t{hap2_tname_possible_suffixes}")
    
    # return hap2_qname_possible_suffixes[0]
    return hap2_qname_possible_suffixes, hap2_tname_possible_suffixes

def find_best_alignment(longest_alignments: Dict[Tuple[str, str], Alignment]) -> Tuple[Tuple[str, str], Alignment]:
    """
    Find the best alignment across all longest alignments based on maximum number of matches.
    
    Returns: ((query_path, ref_path), best_alignment) tuple
    """
    if not longest_alignments:
        raise ValueError("No alignments provided")
    
    # our best alignment is the one with the earliest query start position (we consider all alignments that start within 500 bases of the query start position),
    # and to resolve ties, we select the one with the maximum number of matches
    qstart_sorted_longest_alignments = []
    for (q_path, t_path), aln in longest_alignments.items():
        qstart_sorted_longest_alignments.append((q_path, t_path, aln))
    if len(qstart_sorted_longest_alignments) == 0:
        raise ValueError("No alignments provided")
    qstart_sorted_longest_alignments = sorted(qstart_sorted_longest_alignments, key=lambda x: x[2].q_start)
    candidate_alignments = [qstart_sorted_longest_alignments[0]]
    for alignment in qstart_sorted_longest_alignments[1:]:
        if (alignment[2].q_start - candidate_alignments[0][2].q_start) <= 500:
            candidate_alignments.append(alignment)
    best_alignment = max(candidate_alignments, key=lambda x: x[2].num_matches)
    best_key = (best_alignment[0], best_alignment[1])
    best_aln = best_alignment[2]
    return best_key, best_aln


def find_best_alignment_per_query_ref_combination(alignments: List[Alignment]) -> Dict[Tuple[str, str], Alignment]:
    """
    For each query-path:ref-path combination, find the longest alignment.
    
    Uses alignment block length as the metric for "longest".
    If multiple alignments have the same length, keeps the one with more matches.
    
    Returns: Dictionary mapping (query_path, ref_path) -> longest Alignment
    """
    # Group alignments by (query_path, ref_path) combination
    combinations = defaultdict(list)
    for aln in alignments:
        key = (aln.q_name, aln.t_name)
        combinations[key].append(aln)
    
    longest_alignments = {}
    for (q_path, t_path), alns in combinations.items():
        # Sort by query start position in the alignment (ascending), then by num_matches (descending)
        qstart_sorted_alignments = sorted(alns, key=lambda a: a.q_start)
        # Select all alignments that start within 500 bases of the query start position
        # Among these, select the longest alignment by num_matches
        if len(qstart_sorted_alignments) > 0:
            candidate_alignments = [qstart_sorted_alignments[0]]
            for alignment in qstart_sorted_alignments[1:]:
                if (alignment.q_start - qstart_sorted_alignments[0].q_start) <= 500:
                    candidate_alignments.append(alignment)
            best_alignment = max(candidate_alignments, key=lambda a: a.num_matches)
            longest_alignments[(q_path, t_path)] = best_alignment
    
    return longest_alignments

def find_best_alignment_for_hap2(hap2_possible_query_suffixes, hap2_possible_ref_suffixes, longest_alignments: Dict[Tuple[str, str], Alignment]):
    """
    Find the best alignment for hap2.
    """
    # print(f"longest_alignments: \n\t{longest_alignments.keys()}")
    hap2_alignments = {}
    hap2_primary_alignments = {}

    for hap2_query_suffix in hap2_possible_query_suffixes:
        for hap2_ref_suffix in hap2_possible_ref_suffixes:
            for (q_path, t_path), aln in longest_alignments.items():
                # only consider primary query-path:ref-path combinations
                # an alignment is primary if the 16th (0-indexed) value in it's aln.line is "tp:A:P", and not "tp:A:S".
                # we need to split the aln.line by '\t' and then check the 16th value.
                aln_line_split = aln.line.split('\t')
                primary_alignment_str = aln_line_split[16]
                is_primary_alignment = (primary_alignment_str == "tp:A:P")
                if is_primary_alignment:
                    if q_path.split(':')[1] == hap2_query_suffix and t_path.split(':')[1] == hap2_ref_suffix:
                        hap2_primary_alignments[(q_path, t_path)] = aln
                        break

    # Construct all possible qname, tname pairs and
    # Extract alignments with these qname, tname pairs
    for hap2_query_suffix in hap2_possible_query_suffixes:
        for hap2_ref_suffix in hap2_possible_ref_suffixes:
            for (q_path, t_path), aln in longest_alignments.items():
                if q_path.split(':')[1] == hap2_query_suffix and t_path.split(':')[1] == hap2_ref_suffix:
                    hap2_alignments[(q_path, t_path)] = aln
                    break

    # Find the best alignment for hap2
    # print(f"  Hap2 primary alignments: \n\t{hap2_primary_alignments}")
    # print(f"  Hap2 alignments: \n\t{hap2_alignments}")
    best_key_hap2, best_aln_hap2 = None, None
    if hap2_primary_alignments:
        best_key_hap2, best_aln_hap2 = find_best_alignment(hap2_primary_alignments)
    if ((best_key_hap2 is None) or (best_aln_hap2 is None)):
        best_key_hap2, best_aln_hap2 = find_best_alignment(hap2_alignments)
    # print(f"  Best alignment for hap2: {best_key_hap2}, {best_aln_hap2}")
    return best_key_hap2, best_aln_hap2


def find_overlap_end_in_query_path(query_path:str, query_start_in_aln:int, query_end_in_aln:int, query_contig_graph:Dict[str, Contig]):
    """
    Find the overlap end in the query path.
    We need to find the last contig in the query path that is part of the alignment.
    And then, we need to find the end offset of the alignment in that contig.
    Return the end offset.
    """

    query_path_suffix = query_path.split(':')[1]
    query_path_contig_names = query_path_suffix.split('_')
    query_path_end_contig_name = None
    query_path_end_offset_in_end_contig = -1
    remaining_query_len = query_end_in_aln
    is_first_contig = True
    for contig_idx, contig_name in enumerate(query_path_contig_names):
        curr_contig_overlap_len = 0
        if is_first_contig:
            is_first_contig = False
        else:
            # skip over the overlap length with the previous contig
            last_contig_name = query_path_contig_names[contig_idx - 1]
            for link in query_contig_graph.get(last_contig_name, []).links:
                if link[0] == contig_name:
                    curr_contig_overlap_len = link[2]
                    break
        curr_contig_len = query_contig_graph[contig_name].length
        if remaining_query_len <= curr_contig_len - curr_contig_overlap_len:
            query_path_end_contig_name = contig_name
            query_path_end_offset_in_end_contig = remaining_query_len + curr_contig_overlap_len
            break
        remaining_query_len -= (curr_contig_len - curr_contig_overlap_len)
    
    if query_path_end_contig_name is None:
        raise ValueError(f"Query path end contig not found in query path: {query_path}. This means that the query path is not long enough to contain the alignment. Investigate!!!")

    return query_path_end_contig_name, query_path_end_offset_in_end_contig


def dfs_over_contig_graph(contig_graph: Dict[str, Contig], current_contig_name: str, visited_contig_list: List[str]):

    if current_contig_name not in contig_graph:
        return
    visited_contig_list.append(current_contig_name)
    for link in contig_graph[current_contig_name].links:
        if (link[0] not in visited_contig_list) and (link[0] in contig_graph):
            dfs_over_contig_graph(contig_graph, link[0], visited_contig_list)


def construct_combined_gfa(target_chunk_id, query_chunk_id, target_gfa_file_path, query_gfa_file_path, stitched_gfa_file_path, target_contig_graph, query_contig_graph, query_path_end_contigs_boundaries, all_query_contig_names, target_end_contigs, hap1_aln_record, hap2_aln_record):
    """Construct combined GFA and return case type string."""
    case_type = None
    # NOTE(copied from function call location): 
    # ASSUMPTIONS:
    # 1. This function assumes that the overlap end contig in query paths will be the last contig in the query path itself. (this should be the case if the chunk-boundary region is good) 
    #   open query_gfa_path and target_gfa_path files in read mode, and stitched_gfa_path file in write mode
    # 2. Currently, we don't handle the scenario where the overlap end offset in the overlap end contig leaves less room in the end contig than the length of link overlap with it's next contig. 
    #   If that scenario is hit, there could be erroneous results without warning/exceptions
    # 3. Currently, we only correctly handle scenarios where the to_orientation and from_orientation of the links are both '+'.
    #   If that isn't the case, the stitching might or might not be correct, based on which scenario is hit.

    # IMPLEMENTATION NOTES:
    # 1. Each contig is prepended with its chunk id. (for example, if the chunk id is 52, the contig name will be 52_contig_name)
    # This is done to ensure that the contig names are unique across the stitched GFA file.
    with open(query_gfa_file_path, 'r') as query_gfa_file:
        with open(target_gfa_file_path, 'r') as target_gfa_file:
            with open(stitched_gfa_file_path, 'w') as stitched_gfa_file:
                # read the query_gfa_path file line by line
                for line in target_gfa_file:
                    # write the 'S' lines to the stitched_gfa_path file that don't contain any of the target_end_contigs
                    if line.startswith('S'):
                        line_split = line.strip().split('\t')
                        contig_name = line_split[1]
                        new_contig_name = f"{target_chunk_id}#{contig_name}"
                        line_split[1] = new_contig_name
                        if contig_name not in target_end_contigs:
                            line_remade = '\t'.join(line_split)
                            stitched_gfa_file.write(line_remade + "\n")
                
                print(f"\n\nall_query_contig_names: {all_query_contig_names}\n\n")
                for line in query_gfa_file:
                    # write the 'S' lines to the stitched_gfa_path file that don't contain any of the all_query_contig_names (as they are part of the overlap, and will be cut off. Only the end contig(s) in the query path(s) will be kept, which will be joined to the target overlap end(s) and added later)
                    if line.startswith('S'):
                        line_split = line.strip().split('\t')
                        contig_name = line_split[1]
                        new_contig_name = f"{query_chunk_id}#{contig_name}"
                        line_split[1] = new_contig_name
                        if contig_name not in all_query_contig_names:
                            line_remade = '\t'.join(line_split)
                            stitched_gfa_file.write(line_remade + "\n")
                
                # Now, we concatenate the end contigs of the query path(s) to the target overlap end(s). Note that there are 4 possible scenarios:
                # 1. Query chunk has diploid case, while target chunk is haploid. In this case, we can't concatenate, rather, we would need to add link lines connecting target end contig to the query end contigs.
                # 2. Both query and target chunks are haploid. In this case, we can concatenate the end contigs of the query path(s) to the target overlap end(s).
                # 3. Both query and target chunks are diploid. In this case, we can concatenate the end contigs of the query path(s) to the target overlap end(s). We will need to find the target end contig that matches the query end contig first.
                # 4. Query chunk has haploid case, while target chunk is diploid. In this case, we can't concatenate, rather, we would need to add link lines connecting target end contigs to the query end contig.

                # TODO: Implement scenarios 1. Rest are implemented.
                if (len(query_path_end_contigs_boundaries) == 1) or (len(query_path_end_contigs_boundaries) == 1):
                    raise ValueError("End contig names list cannot be singleton. There must be 2 copies, corresponding to best aln and best aln for hap2, even in haploid case (where they will be duplicates). Investigate!!!")
                if ((query_path_end_contigs_boundaries[0][0] == query_path_end_contigs_boundaries[1][0]) and (target_end_contigs[0] == target_end_contigs[1])):
                    # Scenario 2: Haploid query haploid target
                    case_type = "haploid-haploid"
                    print(f"  Haploid query haploid target. Concatenating end contigs of query path(s) to the target overlap end(s).")
                    query_end_contig_target_end_contig_mapping = [(query_path_end_contigs_boundaries[0][0], query_path_end_contigs_boundaries[0][1], target_end_contigs[0])]
                    
                    # now we concatenate the target end contig and query end contig
                    stitched_contig_records = [] # [(concatenated_contig_name, concatenated_contig_sequence, concatenated_contig_length)]
                    for query_end_contig_name, query_end_contig_offset_in_end_contig, target_end_contig_name in query_end_contig_target_end_contig_mapping:
                        concatenated_contig_name_with_chunk_ids = f"{target_chunk_id}#{target_end_contig_name}+{query_chunk_id}#{query_end_contig_name}"
                        concatenated_contig_name = f"{target_end_contig_name}+{query_end_contig_name}"
                        concatenated_contig_sequence = target_contig_graph[target_end_contig_name].sequence + query_contig_graph[query_end_contig_name].sequence[query_end_contig_offset_in_end_contig:]
                        concatenated_contig_length = len(concatenated_contig_sequence)
                        stitched_gfa_file.write(f"S\t{concatenated_contig_name_with_chunk_ids}\t{concatenated_contig_sequence}\tLN:i:{concatenated_contig_length}\n")
                        stitched_contig_records.append((concatenated_contig_name, concatenated_contig_sequence, concatenated_contig_length))
                    
                    # now we add link lines. iterate over query_gfa_file line by line (and similarly for target_gfa_file) and add link lines for the stitched contig records.
                    query_gfa_file.seek(0)
                    target_gfa_file.seek(0)

                    old_lines_to_replace_by_stitched_lines_in_query = []
                    for line in query_gfa_file:
                        if line.startswith('L'):
                            line_split = line.strip().split('\t')
                            from_contig_name = line_split[1]
                            from_contig_name_with_chunk_id = query_chunk_id + '#' + from_contig_name
                            line_split[1] = from_contig_name_with_chunk_id
                            to_contig_name = line_split[3]
                            to_contig_name_with_chunk_id = query_chunk_id + '#' + to_contig_name
                            line_split[3] = to_contig_name_with_chunk_id
                            
                            if (from_contig_name not in all_query_contig_names) and (to_contig_name not in all_query_contig_names):
                                line_remade = '\t'.join(line_split)
                                # print(f"  Writing line: {line_remade}")
                                stitched_gfa_file.write(line_remade + "\n")
                            elif from_contig_name in [query_path_end_contigs_boundaries[0][0], query_path_end_contigs_boundaries[1][0]]:
                                old_lines_to_replace_by_stitched_lines_in_query.append(line)
                                # print(f"NOT writing line right now: {line.strip()}. Will decorate it with chunk ids later when writing it.")
                    
                    old_lines_to_replace_by_stitched_lines_in_target = []
                    for line in target_gfa_file:
                        if line.startswith('L'):
                            line_split = line.strip().split('\t')
                            from_contig_name = line_split[1]
                            from_contig_name_with_chunk_id = target_chunk_id + '#' + from_contig_name
                            line_split[1] = from_contig_name_with_chunk_id
                            to_contig_name = line_split[3]
                            to_contig_name_with_chunk_id = target_chunk_id + '#' + to_contig_name
                            line_split[3] = to_contig_name_with_chunk_id
                            if (from_contig_name not in target_end_contigs) and (to_contig_name not in target_end_contigs):
                                line_remade = '\t'.join(line_split)
                                # print(f"  Writing line: {line_remade}")
                                stitched_gfa_file.write(line_remade + "\n")
                            elif to_contig_name in target_end_contigs:
                                # print(f"NOT writing line right now: {line.strip()}. Will decorate it with chunk ids later when writing it.")
                                old_lines_to_replace_by_stitched_lines_in_target.append(line)
                    
                    # now we add the stitched contig records to the stitched_gfa_file
                    # first decorating the stitched contig name with chunk ids
                        target_portion_of_stitched_contig_name, query_portion_of_stitched_contig_name = stitched_contig_records[0][0].split('+')
                        stitched_contig_name_with_chunk_ids = f"{target_chunk_id}#{target_portion_of_stitched_contig_name}+{query_chunk_id}#{query_portion_of_stitched_contig_name}"

                    for old_line in old_lines_to_replace_by_stitched_lines_in_query:
                        old_line_split = old_line.strip().split('\t')
                        old_to_contig_name = old_line_split[3]
                        old_to_contig_name_with_chunk_id = query_chunk_id + '#' + old_to_contig_name
                        # skipping the from_contig_name part as we use the stitched_contig_name here
                        old_overlap_len = old_line_split[5]
                        old_from_orientation = old_line_split[2]
                        old_to_orientation = old_line_split[4]
                        stitched_gfa_file.write(f"L\t{stitched_contig_name_with_chunk_ids}\t{old_from_orientation}\t{old_to_contig_name_with_chunk_id}\t{old_to_orientation}\t{old_overlap_len}\n")
                        
                    for old_line in old_lines_to_replace_by_stitched_lines_in_target:
                        old_line_split = old_line.strip().split('\t')
                        old_from_contig_name = old_line_split[1]
                        old_from_contig_name_with_chunk_id = target_chunk_id + '#' + old_from_contig_name
                        old_overlap_len = old_line_split[5]
                        old_from_orientation = old_line_split[2]
                        old_to_orientation = old_line_split[4]
                        stitched_gfa_file.write(f"L\t{old_from_contig_name_with_chunk_id}\t{old_from_orientation}\t{stitched_contig_name_with_chunk_ids}\t{old_to_orientation}\t{old_overlap_len}\n")
                
                elif query_path_end_contigs_boundaries[0][0] == query_path_end_contigs_boundaries[1][0]:
                    case_type = "haploid-diploid"
                    raise NotImplementedError("This 'haploid query:diploid target' scenario is not implemented yet.")
                
                elif ((target_end_contigs[0] == target_end_contigs[1]) and (query_path_end_contigs_boundaries[0][0] != query_path_end_contigs_boundaries[1][0])):
                    
                    #Scenario 4: Diploid query haploid target.
                    case_type = "diploid-haploid"
                    print(f"  Overlap end-contig(s) in query are diploid, while target end-contig is haploid. Linking target end-contig to both the query overlap end-contig(s), after trimming the overlap end-contig(s) in the query")
                    # first we find which query end contig maps to which target end contig
                    print(f"  Target end contig(s): {target_end_contigs}")
                    query_end_contig_target_end_contig_mapping = []  # [(query_end_contig_name, query_end_contig_offset_in_end_contig, target_end_contig_name)]
                    hap1_query_end_contig_idx = 1
                    hap1_query_contigs = hap1_aln_record.q_name.split(':')[1].split('_')
                    if query_path_end_contigs_boundaries[0][0] in hap1_query_contigs:
                        hap1_query_end_contig_idx = 0
                    query_end_contig_target_end_contig_mapping.append((query_path_end_contigs_boundaries[0][0], query_path_end_contigs_boundaries[0][1], target_end_contigs[0]))
                    query_end_contig_target_end_contig_mapping.append((query_path_end_contigs_boundaries[1][0], query_path_end_contigs_boundaries[1][1], target_end_contigs[1]))

                    # since we're not stitching contigs in this case, we add the overlap-end contigs from query & target end contigs back to the stitched_gfa_file
                    target_end_contig_name_with_chunk_id = target_chunk_id + '#' + target_end_contigs[0]
                    stitched_gfa_file.write(f"S\t{target_end_contig_name_with_chunk_id}\t{target_contig_graph[target_end_contigs[0]].sequence}\tLN:i:{target_contig_graph[target_end_contigs[0]].length}\n")
                    stitched_contig_records = [] # [(concatenated_contig_name, concatenated_contig_sequence, concatenated_contig_length)]
                    
                    for query_end_contig_name, query_end_contig_offset_in_end_contig, _ in query_end_contig_target_end_contig_mapping:
                        query_contig_name = f"{query_end_contig_name}-tr"
                        query_contig_name_with_chunk_id = query_chunk_id + '#' + query_contig_name
                        query_contig_sequence = query_contig_graph[query_end_contig_name].sequence[query_end_contig_offset_in_end_contig:]
                        query_contig_length = len(query_contig_sequence)
                        stitched_gfa_file.write(f"S\t{query_contig_name_with_chunk_id}\t{query_contig_sequence}\tLN:i:{query_contig_length}\n")
                        stitched_contig_records.append((query_contig_name, query_contig_sequence, query_contig_length))
                    
                    # now we add link lines. iterate over query_gfa_file line by line (and similarly for target_gfa_file) and add link lines for the stitched contig records.
                    query_gfa_file.seek(0)
                    target_gfa_file.seek(0)

                    # now we link the query overlap end contigs to the target end contig
                    old_lines_to_replace_by_trimmed_lines_in_query = []
                    for line in query_gfa_file:
                        if line.startswith('L'):
                            line_split = line.strip().split('\t')
                            from_contig_name = line_split[1]
                            from_contig_name_with_chunk_id = query_chunk_id + '#' + from_contig_name
                            line_split[1] = from_contig_name_with_chunk_id
                            to_contig_name = line_split[3]
                            to_contig_name_with_chunk_id = query_chunk_id + '#' + to_contig_name
                            line_split[3] = to_contig_name_with_chunk_id
                            if (from_contig_name not in all_query_contig_names) and (to_contig_name not in all_query_contig_names):
                                line_remade = '\t'.join(line_split)
                                stitched_gfa_file.write(line_remade + "\n")
                            elif from_contig_name in [query_path_end_contigs_boundaries[0][0], query_path_end_contigs_boundaries[1][0]]:
                                old_lines_to_replace_by_trimmed_lines_in_query.append(line)
                    
                    for line in target_gfa_file:
                        if line.startswith('L'):
                            line_split = line.strip().split('\t')
                            from_contig_name = line_split[1]
                            from_contig_name_with_chunk_id = target_chunk_id + '#' + from_contig_name
                            line_split[1] = from_contig_name_with_chunk_id
                            to_contig_name = line_split[3]
                            to_contig_name_with_chunk_id = target_chunk_id + '#' + to_contig_name
                            line_split[3] = to_contig_name_with_chunk_id
                            if (to_contig_name in target_end_contigs) or (from_contig_name in target_end_contigs):
                                line_remade = '\t'.join(line_split)
                                print(f"Skipping writing line: {line_remade}")
                            if (from_contig_name not in target_end_contigs):    # Note that in this scenario, we aren't changing the target end contig itself, so we can also keep the original links where target end contig is a sink.
                                line_remade = '\t'.join(line_split)
                                stitched_gfa_file.write(line_remade + "\n")
                    
                    # adding (to the stitched_gfa_file) the links originating from the trimmed query overlap end contig
                    for old_line in old_lines_to_replace_by_trimmed_lines_in_query:
                        old_line = old_line.strip().split('\t')
                        old_to_contig_name = old_line[3]
                        old_to_contig_name_with_chunk_id = query_chunk_id + '#' + old_to_contig_name
                        old_from_contig_name = old_line[1]
                        old_overlap_len_str = old_line[5]
                        old_overlap_len = int(old_overlap_len_str[:-1])
                        old_from_orientation = old_line[2]
                        old_to_orientation = old_line[4]
                        idx_in_query_end_contig_mapping = 0 if (old_from_contig_name == query_end_contig_target_end_contig_mapping[0][0]) else 1
                        stitched_contig_name_with_chunk_ids = f"{query_chunk_id}#{stitched_contig_records[idx_in_query_end_contig_mapping][0]}"

                        if stitched_contig_records[idx_in_query_end_contig_mapping][2] < old_overlap_len:
                            raise ValueError(f"Trimmed query ovelap end contig is shorter than old link overlap length ({stitched_contig_records[idx_in_query_end_contig_mapping][2]} vs {old_overlap_len}). This case isn't handled yet.")
                        stitched_gfa_file.write(f"L\t{stitched_contig_name_with_chunk_ids}\t{old_from_orientation}\t{old_to_contig_name_with_chunk_id}\t{old_to_orientation}\t{old_overlap_len_str}\n")

                    # adding (to the stitched gfa file) the links originating from target end contig and sinking at query overlap end contig(s)
                    for query_overlap_end_contig_name, _, _ in stitched_contig_records:
                        query_overlap_end_contig_name_with_chunk_id = query_chunk_id + '#' + query_overlap_end_contig_name
                        target_end_contig_name_with_chunk_id = target_chunk_id + '#' + target_end_contigs[0]
                        stitched_gfa_file.write(f"L\t{target_end_contig_name_with_chunk_id}\t+\t{query_overlap_end_contig_name_with_chunk_id}\t+\t0M\n") 
                
                elif ((target_end_contigs[0] != target_end_contigs[1]) and (query_path_end_contigs_boundaries[0][0] != query_path_end_contigs_boundaries[1][0])):
                    #Scenario 3: Both query and target chunks are diploid.
                    case_type = "diploid-diploid"
                    print(f"  Both query and target chunks are diploid. Concatenating end contigs of query path(s) to the target overlap end(s).")
                    # first we find which query end contig maps to which target end contig
                    query_end_contig_target_end_contig_mapping = []  # [(query_end_contig_name, query_end_contig_offset_in_end_contig, target_end_contig_name)]
                    hap1_query_end_contig_idx = 1
                    hap1_query_contigs = hap1_aln_record.q_name.split(':')[1].split('_')
                    if query_path_end_contigs_boundaries[0][0] in hap1_query_contigs:
                        hap1_query_end_contig_idx = 0
                    query_end_contig_target_end_contig_mapping.append((query_path_end_contigs_boundaries[hap1_query_end_contig_idx][0], query_path_end_contigs_boundaries[hap1_query_end_contig_idx][1], hap1_aln_record.t_name.split(':')[1].split('_')[-1]))

                    query_end_contig_target_end_contig_mapping.append((query_path_end_contigs_boundaries[1 - hap1_query_end_contig_idx][0], query_path_end_contigs_boundaries[1 - hap1_query_end_contig_idx][1], hap2_aln_record.t_name.split(':')[1].split('_')[-1]))

                    # now we concatenate the respective end contigs
                    stitched_contig_records = [] # [(concatenated_contig_name, concatenated_contig_sequence, concatenated_contig_length)]
                    for query_end_contig_name, query_end_contig_offset_in_end_contig, target_end_contig_name in query_end_contig_target_end_contig_mapping:
                        concatenated_contig_name_with_chunk_ids = f"{target_chunk_id}#{target_end_contig_name}+{query_chunk_id}#{query_end_contig_name}"
                        concatenated_contig_name = f"{target_end_contig_name}+{query_end_contig_name}"
                        concatenated_contig_sequence = target_contig_graph[target_end_contig_name].sequence + query_contig_graph[query_end_contig_name].sequence[query_end_contig_offset_in_end_contig:]
                        concatenated_contig_length = len(concatenated_contig_sequence)
                        stitched_gfa_file.write(f"S\t{concatenated_contig_name_with_chunk_ids}\t{concatenated_contig_sequence}\tLN:i:{concatenated_contig_length}\n")
                        stitched_contig_records.append((concatenated_contig_name, concatenated_contig_sequence, concatenated_contig_length))
                    
                    # now we add link lines. iterate over query_gfa_file line by line (and similarly for target_gfa_file) and add link lines for the stitched contig records.
                    query_gfa_file.seek(0)
                    target_gfa_file.seek(0)

                    old_lines_to_replace_by_stitched_lines_in_query = []
                    for line in query_gfa_file:
                        if line.startswith('L'):
                            line_split = line.strip().split('\t')
                            from_contig_name = line_split[1]
                            from_contig_name_with_chunk_id = query_chunk_id + '#' + from_contig_name
                            line_split[1] = from_contig_name_with_chunk_id
                            to_contig_name = line_split[3]
                            to_contig_name_with_chunk_id = query_chunk_id + '#' + to_contig_name
                            line_split[3] = to_contig_name_with_chunk_id
                            if (from_contig_name not in all_query_contig_names) and (to_contig_name not in all_query_contig_names):
                                line_remade = '\t'.join(line_split)
                                stitched_gfa_file.write(line_remade + "\n")
                            elif from_contig_name in [query_path_end_contigs_boundaries[0][0], query_path_end_contigs_boundaries[1][0]]:
                                old_lines_to_replace_by_stitched_lines_in_query.append(line)
                    
                    old_lines_to_replace_by_stitched_lines_in_target = []
                    for line in target_gfa_file:
                        if line.startswith('L'):
                            line_split = line.strip().split('\t')
                            from_contig_name = line_split[1]
                            from_contig_name_with_chunk_id = target_chunk_id + '#' + from_contig_name
                            line_split[1] = from_contig_name_with_chunk_id
                            to_contig_name = line_split[3]
                            to_contig_name_with_chunk_id = target_chunk_id + '#' + to_contig_name
                            line_split[3] = to_contig_name_with_chunk_id
                            if (from_contig_name not in target_end_contigs) and (to_contig_name not in target_end_contigs):
                                line_remade = '\t'.join(line_split)
                                stitched_gfa_file.write(line_remade + "\n")
                            elif to_contig_name in target_end_contigs:
                                old_lines_to_replace_by_stitched_lines_in_target.append(line)
                    
                    # now we add the stitched contig records to the stitched_gfa_file
                    for old_line in old_lines_to_replace_by_stitched_lines_in_query:
                        old_line_split = old_line.strip().split('\t')
                        old_to_contig_name = old_line_split[3]
                        old_to_contig_name_with_chunk_id = query_chunk_id + '#' + old_to_contig_name
                        old_from_contig_name = old_line_split[1]
                        old_overlap_len = old_line_split[5]
                        old_from_orientation = old_line_split[2]
                        old_to_orientation = old_line_split[4]
                        idx_in_query_end_contig_mapping = 0 if (old_from_contig_name == query_end_contig_target_end_contig_mapping[0][0]) else 1
                        stitched_contig_name = stitched_contig_records[idx_in_query_end_contig_mapping][0]
                        target_portion_of_stitched_contig_name, query_portion_of_stitched_contig_name = stitched_contig_name.split('+')
                        stitched_contig_name_with_chunk_ids = f"{target_chunk_id}#{target_portion_of_stitched_contig_name}+{query_chunk_id}#{query_portion_of_stitched_contig_name}"
                        stitched_gfa_file.write(f"L\t{stitched_contig_name_with_chunk_ids}\t{old_from_orientation}\t{old_to_contig_name_with_chunk_id}\t{old_to_orientation}\t{old_overlap_len}\n")

                    for old_line in old_lines_to_replace_by_stitched_lines_in_target:
                        old_line_split = old_line.strip().split('\t')
                        old_from_contig_name = old_line_split[1]
                        old_from_contig_name_with_chunk_id = target_chunk_id + '#' + old_from_contig_name
                        old_to_contig_name = old_line_split[3]
                        old_overlap_len = old_line_split[5]
                        old_from_orientation = old_line_split[2]
                        old_to_orientation = old_line_split[4]
                        idx_in_target_end_contig_mapping = 0 if (old_to_contig_name == query_end_contig_target_end_contig_mapping[0][2]) else 1
                        stitched_contig_name = stitched_contig_records[idx_in_target_end_contig_mapping][0]
                        target_portion_of_stitched_contig_name, query_portion_of_stitched_contig_name = stitched_contig_name.split('+')
                        stitched_contig_name_with_chunk_ids = f"{target_chunk_id}#{target_portion_of_stitched_contig_name}+{query_chunk_id}#{query_portion_of_stitched_contig_name}"
                        stitched_gfa_file.write(f"L\t{old_from_contig_name_with_chunk_id}\t{old_from_orientation}\t{stitched_contig_name_with_chunk_ids}\t{old_to_orientation}\t{old_overlap_len}\n")
    
    return case_type

def print_summary(longest_alignments: Dict[Tuple[str, str], Alignment]):
    """Print summary statistics."""
    print(f"\nSummary:")
    print(f"  Total query-path:ref-path combinations: {len(longest_alignments)}")
    
    # Count unique query paths and ref paths
    query_paths = set(q for q, _ in longest_alignments.keys())
    ref_paths = set(t for _, t in longest_alignments.keys())
    
    print(f"  Unique query paths: {len(query_paths)}")
    print(f"  Unique ref paths: {len(ref_paths)}")
    
    # Show alignment lengths
    if longest_alignments:
        lengths = [aln.aln_block_len for aln in longest_alignments.values()]
        print(f"  Alignment lengths - Min: {min(lengths)}, Max: {max(lengths)}")

        matches = [aln.num_matches for aln in longest_alignments.values()]
        print(f"  Matches - Min: {min(matches)}, Max: {max(matches)}")



def main():
    parser = argparse.ArgumentParser(
        description='Find longest alignment for each query-path:ref-path combination from PAF file.'
    )
    parser.add_argument('-p', '--paf', required=True, help='Input PAF file with alignments')
    parser.add_argument('-g0', '--gfa0', required=True, help='Input GFA file with Chunk-0')
    parser.add_argument('-g1', '--gfa1', required=True, help='Input GFA file with Chunk-1')
    parser.add_argument('-c', '--chunk_pairs', required=True, help='Underscore delimited target and query chunk ids. Example: 52_53')
    parser.add_argument('-o', '--output', required=True, help='Output filtered PAF file with high scoring alignments only')
    parser.add_argument('-s', '--stitched_gfa', required=True, help='Output stitched GFA file with query and target chunks combined')
    parser.add_argument('-t', '--stats-output', required=True, help='Output TSV file with statistics')
    args = parser.parse_args()
    
    target_chunk_id, query_chunk_id = args.chunk_pairs.split('_')
    _, _ = int(target_chunk_id), int(query_chunk_id)    # only checking if the chunk ids are valid integers (else ValueError will be raised)
    
    print(f"Parsing PAF file: {args.paf}")
    alignments = parse_paf(args.paf)
    print(f"  Found {len(alignments)} total alignments")
    
    if not alignments:
        print("Error: No alignments found in PAF file")
        return 1
    
    print(f"\n=============================================================================")
    print(f"Step 1: Finding longest alignment for each query-path:ref-path combination...")
    print(f"===============================================================================")
    longest_alignments = find_best_alignment_per_query_ref_combination(alignments)
    
    print_summary(longest_alignments)
    
    print(f"\n====================================================================")
    print(f"Step 2: Finding best alignment based on maximum number of matches...")
    print(f"======================================================================")
    best_key, best_aln = find_best_alignment(longest_alignments)
    
    q_path, t_path = best_key
    print(f"\nBest alignment (Hap1 path):")
    print(f"  Query path: {q_path}")
    print(f"  Ref path: {t_path}")
    print(f"  CIGAR: {best_aln.cs_line}")
    print(f"  Alignment block length: {best_aln.aln_block_len:,} bp")
    print(f"  Matches: {best_aln.num_matches:,}")
    print(f"  Query coordinates: {best_aln.q_start:,}-{best_aln.q_end:,} (length: {best_aln.q_len:,})")
    print(f"  Target coordinates: {best_aln.t_start:,}-{best_aln.t_end:,} (length: {best_aln.t_len:,})")
    print(f"  Strand: {best_aln.strand}")
    
    # Calculate match ratio (identity)
    if best_aln.aln_block_len > 0:
        match_ratio = best_aln.num_matches / best_aln.aln_block_len
        print(f"  Match ratio (identity): {match_ratio:.4f} ({match_ratio*100:.2f}%)")
    
    print(f"\n================================================")
    print(f"Step 3: Finding best alignment for Hap2 path...")
    print(f"==================================================")

    query_connectivity_dict, ref_connectivity_dict, sibling_grouped_query_connectivity_dict, sibling_grouped_ref_connectivity_dict = generate_contig_relationships_dict(alignments)
    print(f"  Query connectivity dictionary: \n\t{query_connectivity_dict}\n")
    print(f"  Ref connectivity dictionary: \n\t{ref_connectivity_dict}\n")
    
    # Construct Path Name for Hap2 path (qname and tname)        
    hap2_qname_possible_suffixes, hap2_tname_possible_suffixes = construct_hap2_paths(sibling_grouped_query_connectivity_dict, sibling_grouped_ref_connectivity_dict, best_aln)

    # Compare hap1 and hap2 paths. 
    # If they both are the same, then we need not pick the best alignment for hap2
    hap1_qname_suffix = best_aln.q_name.split(':')[1]
    hap1_tname_suffix = best_aln.t_name.split(':')[1]

    if (
        (len(hap2_qname_possible_suffixes) ==  1 and hap2_qname_possible_suffixes[0] == hap1_qname_suffix)
        and (len(hap2_tname_possible_suffixes) ==  1 and hap2_tname_possible_suffixes[0] == hap1_tname_suffix)
    ):
        best_key_hap2, best_aln_hap2 = best_key, best_aln
        print(f"  Both query and reference paths are haploid paths")
    else:
        # Pick the best alignment for hap2 query and reference paths
        best_key_hap2, best_aln_hap2 = find_best_alignment_for_hap2(hap2_qname_possible_suffixes, hap2_tname_possible_suffixes, longest_alignments)
        print(f"  Best alignment for Hap2 path:")
        print(f"    Query path: {best_key_hap2[0]}")
        print(f"    Ref path: {best_key_hap2[1]}")
        print(f"    CIGAR: {best_aln_hap2.cs_line}")
        print(f"    Alignment block length: {best_aln_hap2.aln_block_len:,} bp")
        print(f"    Matches: {best_aln_hap2.num_matches:,}")
        print(f"    Query coordinates: {best_aln_hap2.q_start:,}-{best_aln_hap2.q_end:,} (length: {best_aln_hap2.q_len:,})")
        print(f"    Target coordinates: {best_aln_hap2.t_start:,}-{best_aln_hap2.t_end:,} (length: {best_aln_hap2.t_len:,})")
        print(f"    Strand: {best_aln_hap2.strand}")
    
    print(f"\n=============================================================")
    print(f"Step 4: Output filtered PAF with Hap1 and Hap2 paths only...")
    print(f"===============================================================")
    with open(args.paf, 'r') as f:
        with open(args.output, 'w') as f_out:
            for line in f:            
                if (best_aln.line == line or best_aln_hap2.line == line):
                    f_out.write(line)    

    print(f"  Filtered PAF file saved to: {args.output}")


    print(f"\n===================================================")
    print(f"Step 5: Find overlap end(s) in query path(s)...")
    print(f"=====================================================")

    all_query_path_names = [qname.split(':')[1].split('_') for qname in [best_aln.q_name, best_aln_hap2.q_name]]
    all_query_contig_names = list(set([ele for qname in all_query_path_names for ele in qname]))
    all_ref_path_names = [tname.split(':')[1].split('_') for tname in [best_aln.t_name, best_aln_hap2.t_name]]
    all_ref_contig_names = list(set([ele for tname in all_ref_path_names for ele in tname]))
    
    all_query_initial_candidate_contig_names = []
    # iterate over query_connectivity_dict and add the contig names to all_query_initial_candidate_contig_names
    for pos_idx, contig_list in query_connectivity_dict.items():
        all_query_initial_candidate_contig_names.extend(contig_list)
    all_query_initial_candidate_contig_names = list(set(all_query_initial_candidate_contig_names))

    print(f" all_query_contig_names: {all_query_contig_names}")
    gfa0: Dict[str, Contig] = read_gfa_file_for_candidate_contigs(args.gfa0, set(all_ref_contig_names))  # target contigs graph
    gfa1: Dict[str, Contig] = read_gfa_file_for_candidate_contigs(args.gfa1, set(all_query_contig_names))  # query contigs graph
    
    print(f"\n\n\n gfa1: {gfa1}\n\n\n")
    print(f"  GFA target contigs: \n\t{gfa0}\n")
    print(f"  GFA query contigs: \n\t{gfa1}\n")

    # Find the overlap end(s) in the query path(s)
    hap1_query_path_overlap_end_contig_name, hap1_query_path_end_offset_in_end_contig = find_overlap_end_in_query_path(best_aln.q_name, best_aln.q_start, best_aln.q_end, gfa1)
    print(f"  Hap1 overlap end in query path: {hap1_query_path_overlap_end_contig_name} and alignment end offset in end contig: {hap1_query_path_end_offset_in_end_contig}")
    hap1_query_path_contigs = best_aln.q_name.split(':')[1].split('_')
    hap1_query_path_actual_end_contig_name = hap1_query_path_contigs[-1]
    # TODO: Remove this exception raising completely after testing the recent implementation  
    # if hap1_query_path_overlap_end_contig_name != hap1_query_path_actual_end_contig_name:
    #     raise NotImplementedError("Hap1 query path overlap end contig name does not match the actual end contig name. This should not happen if the chunk-boundary region is good.")

    hap2_query_path_overlap_end_contig_name, hap2_query_path_end_offset_in_end_contig = find_overlap_end_in_query_path(best_aln_hap2.q_name, best_aln_hap2.q_start, best_aln_hap2.q_end, gfa1)
    print(f"  Hap2 overlap end in query path: {hap2_query_path_overlap_end_contig_name} and alignment end offset in end contig: {hap2_query_path_end_offset_in_end_contig}")
    hap2_query_path_contigs = best_aln_hap2.q_name.split(':')[1].split('_')
    hap2_query_path_actual_end_contig_name = hap2_query_path_contigs[-1]
    # TODO: Remove this exception raising completely after testing the recent implementation  
    # if hap2_query_path_overlap_end_contig_name != hap2_query_path_actual_end_contig_name:
    #     raise NotImplementedError("Hap2 query path overlap end contig name does not match the actual end contig name. This should not happen if the chunk-boundary region is good.")

    # in case the query overlap contigs are not the end contigs in the query path, we need to not reject them in our stitched GFA file.
    # So, we list the contigs that occur to the right of the query overlap end contigs, and remove them from the rejected query contigs list, i.e., all_query_contig_names.
    # To find these contigs, we first build the candidate contigs graph for the query, and then chase the links originating from query overlap end contigs to the right, and add whatever contigs are encountered.
    gfa1_initial_candidate_contigs: Dict[str, Contig] = read_gfa_file_for_candidate_contigs(args.gfa1, set(all_query_initial_candidate_contig_names))
    # print(f"  GFA1 initial candidate contigs: \n\t{gfa1_initial_candidate_contigs}\n")
    query_contigs_lying_to_right_of_overlap_end_contigs = []
    # performing dfs over query candidate contig graph, starting from the query overlap end contigs, and adding all contigs that are encountered to the query_contigs_lying_to_right_of_overlap_end_contigs list.
    # print(f"  All query initial candidate contig names: {all_query_initial_candidate_contig_names}")
    dfs_over_contig_graph(gfa1_initial_candidate_contigs, hap1_query_path_overlap_end_contig_name, query_contigs_lying_to_right_of_overlap_end_contigs)
    dfs_over_contig_graph(gfa1_initial_candidate_contigs, hap2_query_path_overlap_end_contig_name, query_contigs_lying_to_right_of_overlap_end_contigs)
    # removing the query overlap end contigs from the query_contigs_lying_to_right_of_overlap_end_contigs list, since we don't want to reject them in our stitched GFA file, but they were added during the dfs.
    query_contigs_lying_to_right_of_overlap_end_contigs.remove(hap1_query_path_overlap_end_contig_name)
    query_contigs_lying_to_right_of_overlap_end_contigs.remove(hap2_query_path_overlap_end_contig_name)
    # print(f"  Query contigs lying to the right of overlap end contigs: {query_contigs_lying_to_right_of_overlap_end_contigs}")
    all_query_candidate_contigs_till_overlap_end_contigs = all_query_initial_candidate_contig_names
    for contig in query_contigs_lying_to_right_of_overlap_end_contigs:
        if contig in all_query_candidate_contigs_till_overlap_end_contigs:
            all_query_candidate_contigs_till_overlap_end_contigs.remove(contig)
    # print(f"  All query initial candidate contig names list after removing contigs lying to the right of overlap end contigs: {all_query_initial_candidate_contig_names}")

    # TODO: Like we found overlap end(s) in query path(s), we need to find overlap end(s) in target path(s)
    # Currently, we are assuming that target paths will be fully contained in the overlap, and so we can simply join the end of the end contigs of target path to the overlap end of the query path. 
    
    print(f"\n===================================================")
    print(f"Step 6: Construct combined GFA for query and target chunks...")
    print(f"=====================================================")

    target_end_contigs = [tname.split(':')[1].split('_')[-1] for tname in [best_aln.t_name, best_aln_hap2.t_name]]
    # NOTE: This construct_combined_gfa function assumes that the overlap end contig in query paths will be the last contig in the query path itself. (this should be the case if the chunk-boundary region is good) 

    case_type = construct_combined_gfa(target_chunk_id, query_chunk_id, args.gfa0, args.gfa1, args.stitched_gfa, gfa0, gfa1, [(hap1_query_path_overlap_end_contig_name, hap1_query_path_end_offset_in_end_contig), (hap2_query_path_overlap_end_contig_name, hap2_query_path_end_offset_in_end_contig)], all_query_candidate_contigs_till_overlap_end_contigs, target_end_contigs, best_aln, best_aln_hap2)
    
    # Collect statistics for TSV output
    # 1. Case type - already collected from construct_combined_gfa
    if case_type is None:
        case_type = "unknown"  # Fallback if case_type wasn't set (shouldn't happen)
    
    # 2. Number of subgroups per index in sibling_grouped_query_connectivity_dict
    query_subgroups_per_index = []
    if sibling_grouped_query_connectivity_dict:
        max_idx = max(sibling_grouped_query_connectivity_dict.keys())
        for idx in range(max_idx + 1):
            if idx in sibling_grouped_query_connectivity_dict:
                query_subgroups_per_index.append(len(sibling_grouped_query_connectivity_dict[idx]))
            else:
                query_subgroups_per_index.append(0)
    
    # 3. Number of subgroups per index in sibling_grouped_ref_connectivity_dict
    ref_subgroups_per_index = []
    if sibling_grouped_ref_connectivity_dict:
        max_idx = max(sibling_grouped_ref_connectivity_dict.keys())
        for idx in range(max_idx + 1):
            if idx in sibling_grouped_ref_connectivity_dict:
                ref_subgroups_per_index.append(len(sibling_grouped_ref_connectivity_dict[idx]))
            else:
                ref_subgroups_per_index.append(0)
    
    # 4. Check if hap1 and hap2 overlap end contigs match actual end contigs
    hap1_end_match = (hap1_query_path_overlap_end_contig_name == hap1_query_path_actual_end_contig_name)
    hap2_end_match = (hap2_query_path_overlap_end_contig_name == hap2_query_path_actual_end_contig_name)
    
    # Write statistics to TSV file
    with open(args.stats_output, 'w') as stats_out:
        # Write header
        stats_out.write("case_type\tquery_subgroups_per_index\tref_subgroups_per_index\thap1_hap2_end_matches\n")
        
        # Format query_subgroups_per_index as list
        query_subgroups_str = "[" + ",".join(str(x) for x in query_subgroups_per_index) + "]"
        
        # Format ref_subgroups_per_index as list
        ref_subgroups_str = "[" + ",".join(str(x) for x in ref_subgroups_per_index) + "]"
        
        # Format hap1_hap2_end_matches as tuple
        hap1_hap2_end_matches_str = f"({hap1_end_match},{hap2_end_match})"
        
        # Write data row
        stats_out.write(f"{case_type}\t{query_subgroups_str}\t{ref_subgroups_str}\t{hap1_hap2_end_matches_str}\n")
    
    print(f"\nStatistics written to: {args.stats_output}")

    return 0


if __name__ == "__main__":
    exit(main())

