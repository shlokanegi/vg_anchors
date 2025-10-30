import json
import time
import pickle
import tempfile
import subprocess
import os.path
from sys import stderr, stdout, exit
from collections import defaultdict
import copy
import multiprocessing
from typing import Union 
import assembler.helpers as helpers
from line_profiler import profile as line_profile
# import shasta2

from bdsg.bdsg import PackedGraph
from assembler.anchor import Anchor
from assembler.node import Node
from assembler.config import settings
from assembler.anchor_coverage import AnchorCoverage
from assembler.gtest import GTest
from functools import cmp_to_key

from assembler.read import Read

shared_align_anchor = None


def init_worker_snarl():
    """
    Initializer for the snarl processing multiprocessing pool.
    On Linux with fork, the global shared_align_anchor is already available via copy-on-write.
    """
    pass


def process_each_snarl_chunk_in_worker(chunk_snarl_list: list):
    """
    Worker function to run reliable snarl finding on each snarl chunk
    """
    global shared_align_anchor

    if settings.DEBUG:
        print(f"######### FINDING RELIABLE SNARLS IN CURRENT CHUNK #########")

    t0 = time.time()
    # Create valid anchors within the current chunk
    valid_anchors_in_current_chunk = []

    for snarl_id in chunk_snarl_list:
        anchors = shared_align_anchor.snarl_to_anchors_dictionary[snarl_id]
        for anchor in anchors:
            # TODO: Check memory address of an anchor if it matches the one in the shared memory
            read_info = [
                [read[0], read[1], read[2], read[3]]
                for read in anchor.bp_matched_reads
            ]
            valid_anchors_in_current_chunk.append([anchor, read_info])

    result = shared_align_anchor.find_reliable_snarls(
        valid_anchors=valid_anchors_in_current_chunk, snarl_list=chunk_snarl_list
    )

    if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
        print(f".. Processed {len(chunk_snarl_list)} snarls in {time.time() - t0}s", flush=True, file=stderr)

    return result


class AlignAnchor:

    def __init__(self, threads: int, read_id_map: dict = None) -> None:
        # useful initialization objects
        self.threads = threads
        self.read_id_map = read_id_map
        self.graph = None
        self.snarl_to_anchors_dictionary = defaultdict(list)
        # This dictionary contains all the snarl IDs, i.e. the primary ones as well as the ones made after merging.
        self.snarl_ids_sorted = []
        # This list contains the snarl IDs retained in the current state. For ex - after reliable snarl filtering, this list will contain only the reliable snarls. Similarly, after merging, this list will contain the snarls made after merging and remove the ones that are now merged.
        self.sentinel_to_anchor: dict = dict()
        self.anchor_reads_dict: dict = dict()
        self.next_handle_expand_boundary = None
        # self.anchor_coverage = AnchorCoverage()  # Add coverage tracking
        self.anchor_read_tracking_dict= {} # Add coverage tracking for anchors
        self.independent_anchor_extension_tracking_dict = {}  # For independent anchor extension
        ## For phasing consistency check
        self.reliable_snarls = []
        self.snarl_variant_type_dict = {}
        self.snarl_coverage_dict = {}
        self.snarl_allelic_coverage_dict = {}
        self.linked_snarls_dictionary = {}   # {snarl_id: [linked_snarl_id1, linked_snarl_id2, ...]}
        self.linked_snarls_compatibility_dict = {}   # {snarl_id: {linked_snarl_id1: True/False, linked_snarl_id2: True/False, ...}}
        self.snarl_common_reads_dict = {}
        self.snarl_read_partitions_dict = {} # {primary_snarl: {other_snarl: {"primary": [...], "other": [...]}}}
        # for extended snarls
        self.extended_snarl_coverage_dict = {}
        self.extended_snarl_allelic_coverage_dict = {}

        self.valid_anchors_from_reliable_snarls = []
        self.outputs_for_file = []
        self.runtime_logs = {}
        # self.reads = dict()    # {read_name: read_object}
        self.read_to_snarl_dictionary = {}


    def merge_results(self, result, reads_processed_file_path=None):
        """
        Merges the results from a worker process into the main AlignAnchor instance.
        """

        # Merge the anchor_reads_dict
        for sentinel, anchor_indices in result["anchor_reads_dict"].items():
            for i, reads in anchor_indices.items():
                if reads:
                    self.anchor_reads_dict[sentinel][i].extend(reads)
        
        # Merge the bp_matched_reads back into the main Anchor objects
        for (sentinel, i), reads in result["bp_matched_reads"].items():
            anchor = self.sentinel_to_anchor[sentinel][i]
            anchor.compute_sentinel_bp_length()
            anchor.bp_matched_reads.extend(reads)
        
        if settings.OUTPUT_LOGGING_FILES and reads_processed_file_path is not None:
            if self.read_id_map is None:
                # Append chunk results to the shared TSV; the caller ensures cleanup before first write
                with open(reads_processed_file_path, "a") as f:
                    for read_name, read_data in result["reads_processed"].items():
                        print(f"{read_name}\t{read_data[1]}\t{read_data[2]}\t{read_data[3]}\t{read_data[4]}\t{read_data[5]}\t{read_data[6]}\t{read_data[7]}\t{read_data[8]}\t{read_data[9]}\t{read_data[10]}", file=f)
            else:
                with open(reads_processed_file_path, "a") as f:
                    for read_identifier, read_data in result["reads_processed"].items():
                        line_to_write = f"{read_identifier}\t" + "\t".join(map(str, read_data[1:]))
                        print(line_to_write, file=f)

    def build(self, dict_path: str, packed_graph_path: str) -> None:

        # loading dictionary
        with open(dict_path, 'rb') as in_f:
            self.sentinel_to_anchor = pickle.load(in_f)

        #loading packedgraph
        self.graph = PackedGraph()
        self.graph.deserialize(packed_graph_path)

        # initializing output dictionary
        for sentinel, anchors in self.sentinel_to_anchor.items():
            self.anchor_reads_dict[sentinel] = [[] for _ in range(len(anchors))]


    def readFasta(self, fasta_path: str) -> None:
        self.fasta_path = fasta_path
        self.read_sequences = {}
        with open(fasta_path, "r") as f:
            read_name = None
            for line in f:
                if line.startswith(">"):
                    read_name = line.strip().split()[0][1:]
                    self.read_sequences[read_name] = ""
                elif read_name:
                    self.read_sequences[read_name] += line.strip()
        
        if self.read_id_map:
            self.read_sequences = {self.read_id_map.get(name): seq for name, seq in self.read_sequences.items() if self.read_id_map.get(name) is not None}

    def ingest(self, dictionary: dict, packed_graph_path: str) -> None:
        self.sentinel_to_anchor = dictionary
        self.graph.deserialize(packed_graph_path)

        for sentinel, anchors in self.sentinel_to_anchor.items():
            self.anchor_reads_dict[sentinel] = [[] for _ in range(len(anchors))]


    def _extending_anchors_by_merging(self, snarl_ids_sorted_list_up_to_date, snarl_ids_sorted_list_iterator_idx, current_snarl_id, other_snarl_id, current_snarl_anchors, extend_left, anchors_to_discard, snarl_orientation, merging_round) -> list:
        """
        Attempts to merge anchors from two adjacent snarls. This function is called when we want to combine anchors
        from neighboring snarls to create longer, more robust anchors.

        The function:
        1. Checks if there are enough anchors with sufficient read overlap between the snarls
        2. Creates new merged anchors by combining compatible anchors from both snarls
        3. Updates the snarl IDs and anchor boundaries accordingly
        4. Maintains read coverage information for the merged anchors

        Parameters
        ----------
        snarl_ids_sorted_list_up_to_date : list
            List of all snarl IDs in sorted order
        snarl_ids_sorted_list_iterator_idx : int
            Current position in the snarl IDs list
        current_snarl_id : str
            ID of the current snarl being processed
        other_snarl_id : str
            ID of the adjacent snarl to merge with
        current_snarl_anchors : list
            List of Anchor objects in the current snarl
        extend_left : bool
            True if merging towards left, False if merging towards right
        anchors_to_discard : set
            List to store anchors that should be removed after merging
        snarl_orientation : bool
            Orientation of the snarl (True for forward, False for reverse)
        merging_round : int
            Current round of merging, whether strict-0 (MIN_READS_REQUIRED_FOR_MERGING_R0) or relaxed-1 (MIN_READS_REQUIRED_FOR_MERGING_R1) 
        Returns
        -------
        tuple
            (new_anchors_after_merging, updated_snarl_ids_list_iterator_idx)
            - new_anchors_after_merging: List of newly created merged anchors
            - updated_snarl_ids_list_iterator_idx: Updated position in snarl IDs list
        """

        if settings.DEBUG:
            print("current snarl being extended:- ", current_snarl_id, flush=True, file=stderr)
            print("surrounding snarls:- ", snarl_ids_sorted_list_up_to_date[snarl_ids_sorted_list_iterator_idx-5:min(len(snarl_ids_sorted_list_up_to_date), snarl_ids_sorted_list_iterator_idx+5)], flush=True, file=stderr)
        
        # calculating if after merging, enough anchors (>=2) will remain 
        cnt_anchors_with_sufficient_read_overlap = 0
        # other_anchor = self.snarl_to_anchors_dictionary[other_snarl_id][0]

        # When merging in round=1 (i.e. lower MIN_READS_REQUIRED_FOR_MERGING value, hence less confidence), merge only if current anchor is heterozygous. 
        # Else, there is no point in sacrificing read coverage by such an amount, if it doesn't help in phasing. 
        if merging_round == 1 and len(current_snarl_anchors) < 2:
            return (current_snarl_anchors, snarl_ids_sorted_list_iterator_idx)
        MIN_READS_REQUIRED_FOR_MERGING = settings.MIN_READS_REQUIRED_FOR_MERGING_R0 if merging_round == 0 else settings.MIN_READS_REQUIRED_FOR_MERGING_R1
        
        for other_anchor in self.snarl_to_anchors_dictionary[other_snarl_id]:
            for anchor in current_snarl_anchors:
                common_paths = set(anchor.reference_paths_covered).intersection(set(other_anchor.reference_paths_covered))
                read_ids_current_anchor = [read[settings.READ_POSITION] for read in anchor.bp_matched_reads]
                read_ids_other_anchor = [read[settings.READ_POSITION] for read in other_anchor.bp_matched_reads]
                common_reads_ids = set(read_ids_current_anchor).intersection(set(read_ids_other_anchor))

                # We only consider a common read if it has the same orientation relative to increasing node ids in both anchors.
                read_orientations_relative_to_increasing_node_ids = {} # {read_id: (orientation, [read_start, read_end])}; where orientation = 1 => forward, -1 => reverse
                ORIENTATION_INDEX = 0
                READ_POSITION_INDEX = 1
                READ_START_INDEX = 0
                READ_END_INDEX = 1
                for read in anchor.bp_matched_reads:
                    if read[settings.READ_POSITION] in common_reads_ids:
                        if settings.DEBUG:
                            print(f"DEBUG: In current snarl {current_snarl_id}, read {read[settings.READ_POSITION]} has read info: {read}", flush=True, file=stderr)
                        if read[settings.READ_STRAND] == 0:
                            extra_bps = read[settings.CS_LEFT_AVAIL] if extend_left else read[settings.CS_RIGHT_AVAIL]
                        else:
                            extra_bps = read[settings.CS_RIGHT_AVAIL] if extend_left else read[settings.CS_LEFT_AVAIL]
                        if extra_bps < 1:
                            common_reads_ids.remove(read[settings.READ_POSITION])
                            continue
                        path_orientation = (anchor[0].id < anchor[1].id)
                        read_orientations_relative_to_increasing_node_ids[read[settings.READ_POSITION]] = ((1 if ((read[settings.READ_STRAND] == 0 and path_orientation) or (read[settings.READ_STRAND] == 1 and not path_orientation)) else -1), [read[settings.ANCHOR_START], read[settings.ANCHOR_END]])
                for read in other_anchor.bp_matched_reads:
                    if read[settings.READ_POSITION] in common_reads_ids:
                        if settings.DEBUG:
                            print(f"DEBUG: In other snarl {other_snarl_id}, read {read[settings.READ_POSITION]} has read info: {read}", flush=True, file=stderr)
                        if read[settings.READ_STRAND] == 0:
                            extra_bps = read[settings.CS_RIGHT_AVAIL] if extend_left else read[settings.CS_LEFT_AVAIL]
                        else:
                            extra_bps = read[settings.CS_LEFT_AVAIL] if extend_left else read[settings.CS_RIGHT_AVAIL]
                        if extra_bps < 1:
                            common_reads_ids.remove(read[settings.READ_POSITION])
                            continue
                        path_orientation = (other_anchor[0].id < other_anchor[1].id)
                        current_read_orientation_relative_to_increasing_node_id = 1 if ((read[settings.READ_STRAND] == 0 and path_orientation) or (read[settings.READ_STRAND] == 1 and not path_orientation)) else -1
                        # drop this read if it has conflicting orientations relative to increasing node ids in the two anchors.
                        if (read_orientations_relative_to_increasing_node_ids[read[settings.READ_POSITION]][ORIENTATION_INDEX] != current_read_orientation_relative_to_increasing_node_id):
                            common_reads_ids.remove(read[settings.READ_POSITION])
                            continue
                        sorted_list_of_adjacent_anchor_boundaries_in_read = sorted([read_orientations_relative_to_increasing_node_ids[read[settings.READ_POSITION]][READ_POSITION_INDEX][READ_START_INDEX], read_orientations_relative_to_increasing_node_ids[read[settings.READ_POSITION]][READ_POSITION_INDEX][READ_END_INDEX], read[settings.ANCHOR_START], read[settings.ANCHOR_END]])
                        if sorted_list_of_adjacent_anchor_boundaries_in_read[1] < sorted_list_of_adjacent_anchor_boundaries_in_read[2]:
                            common_reads_ids.remove(read[settings.READ_POSITION])
                            continue

                if settings.DEBUG:
                    print(f"for anchor {anchor!r} and other anchor {other_anchor!r}, common reads between them are {len(common_reads_ids)}", flush=True, file=stderr)                 
                if len(common_paths) > 0 and (len(common_reads_ids) > MIN_READS_REQUIRED_FOR_MERGING):    # meaning we can create an anchor with this combination
                    cnt_anchors_with_sufficient_read_overlap += 1
                    if cnt_anchors_with_sufficient_read_overlap >= 2:
                        break
                    # extendable_anchors_in_current_snarl.append(anchor)
            if cnt_anchors_with_sufficient_read_overlap > 1:
                break
        
        if cnt_anchors_with_sufficient_read_overlap < 2:
            if settings.DEBUG:
                print(f"Failed to merge snarls {current_snarl_id} and {other_snarl_id} because of insufficient anchors after merging", flush=True, file=stderr)
            return (current_snarl_anchors, snarl_ids_sorted_list_iterator_idx)
        
        # meaning we found that after extension, heterozygosity of current snarl will be maintained.
        # extend this side
        new_anchors_after_merging = []
        # loop over all possible combinations of anchor_in_current_snarl x anchor_in_other_snarl
        for anchor in current_snarl_anchors:
            for other_anchor in self.snarl_to_anchors_dictionary[other_snarl_id]:
                common_paths = set(anchor.reference_paths_covered).intersection(set(other_anchor.reference_paths_covered))
                read_ids_current_anchor = [read[settings.READ_POSITION] for read in anchor.bp_matched_reads]
                read_ids_other_anchor = [read[settings.READ_POSITION] for read in other_anchor.bp_matched_reads]
                common_reads_ids = set(read_ids_current_anchor).intersection(set(read_ids_other_anchor))
                
                # We only consider a common read if it has the same orientation relative to increasing node ids in both anchors.
                read_orientations_relative_to_increasing_node_ids = {} # {read_id: (orientation, [read_start, read_end])}; where orientation = 1 => forward, -1 => reverse
                ORIENTATION_INDEX = 0
                READ_POSITION_INDEX = 1
                READ_START_INDEX = 0
                READ_END_INDEX = 1

                for read in anchor.bp_matched_reads:
                    if read[settings.READ_POSITION] in common_reads_ids:
                        if read[settings.READ_STRAND] == 0:
                            extra_bps = read[settings.CS_LEFT_AVAIL] if extend_left else read[settings.CS_RIGHT_AVAIL]
                        else:
                            extra_bps = read[settings.CS_RIGHT_AVAIL] if extend_left else read[settings.CS_LEFT_AVAIL]
                        if extra_bps < 1:
                            common_reads_ids.remove(read[settings.READ_POSITION])
                            if settings.DEBUG:
                                print(f"... read {read[settings.READ_POSITION]} rejected because of possibility of insertion/deletion.", flush=True, file=stderr)
                            continue
                        path_orientation = (anchor[0].id < anchor[1].id)
                        read_orientations_relative_to_increasing_node_ids[read[settings.READ_POSITION]] = ((1 if ((read[settings.READ_STRAND] == 0 and path_orientation) or (read[settings.READ_STRAND] == 1 and not path_orientation)) else -1), [read[settings.ANCHOR_START], read[settings.ANCHOR_END]])

                for read in other_anchor.bp_matched_reads:
                    if read[settings.READ_POSITION] in common_reads_ids:
                        if read[settings.READ_STRAND] == 0:
                            extra_bps = read[settings.CS_RIGHT_AVAIL] if extend_left else read[settings.CS_LEFT_AVAIL]
                        else:
                            extra_bps = read[settings.CS_LEFT_AVAIL] if extend_left else read[settings.CS_RIGHT_AVAIL]
                        if extra_bps < 1:
                            common_reads_ids.remove(read[settings.READ_POSITION])
                            if settings.DEBUG:
                                print(f"... read {read[settings.READ_POSITION]} rejected because of possibility of insertion/deletion.", flush=True, file=stderr)
                            continue
                        path_orientation = (other_anchor[0].id < other_anchor[1].id)
                        current_read_orientation_relative_to_increasing_node_id = 1 if ((read[settings.READ_STRAND] == 0 and path_orientation) or (read[settings.READ_STRAND] == 1 and not path_orientation)) else -1
                        # drop this read if it has conflicting orientations relative to increasing node ids in the two anchors.
                        if (read_orientations_relative_to_increasing_node_ids[read[settings.READ_POSITION]][ORIENTATION_INDEX] != current_read_orientation_relative_to_increasing_node_id):
                            common_reads_ids.remove(read[settings.READ_POSITION])
                            continue
                        sorted_list_of_adjacent_anchor_boundaries_in_read = sorted([read_orientations_relative_to_increasing_node_ids[read[settings.READ_POSITION]][READ_POSITION_INDEX][READ_START_INDEX], read_orientations_relative_to_increasing_node_ids[read[settings.READ_POSITION]][READ_POSITION_INDEX][READ_END_INDEX], read[settings.ANCHOR_START], read[settings.ANCHOR_END]])
                        if sorted_list_of_adjacent_anchor_boundaries_in_read[1] < sorted_list_of_adjacent_anchor_boundaries_in_read[2]:
                            common_reads_ids.remove(read[settings.READ_POSITION])
                            continue


                if len(common_paths) > 0 and (len(common_reads_ids) > MIN_READS_REQUIRED_FOR_MERGING):    # meaning we can create an anchor with this combination
                    if settings.DEBUG:
                        print(f"merging anchors {anchor!r} (current_anchor) and {other_anchor!r} (other_anchor) in", "left extension." if extend_left==True else "right extension.", end=" ", flush=True, file=stderr)
                        if extend_left:
                            print(f"{anchor!r}.bp_occupied_start_node={anchor.bp_occupied_start_node}, other_anchor.bp_occupied_end_node={other_anchor.bp_occupied_end_node}")
                        else:
                            print(f"{anchor!r}.bp_occupied_end_node={anchor.bp_occupied_end_node}, other_anchor.bp_occupied_start_node={other_anchor.bp_occupied_start_node}")
                    new_anchor = copy.deepcopy(anchor)
                    # new_anchor = anchor
                    # find relative orientations of anchor and other anchor
                    new_anchor_orientation = True if (new_anchor[0].id < new_anchor[1].id) else False
                    other_anchor_orientation = True if (other_anchor[0].id < other_anchor[1].id) else False
                    if new_anchor_orientation != other_anchor_orientation:    # flip other anchor here
                        other_anchor.flip_anchor()
                    insert_left = extend_left
                    if new_anchor_orientation != snarl_orientation:
                        insert_left = not insert_left
                    if not new_anchor.merge_anchor(other_anchor, insert_left=insert_left):
                        return (current_snarl_anchors, snarl_ids_sorted_list_iterator_idx)
                    # new snarl id calculation
                    new_snarl_id = (str(other_anchor.snarl_id) + "-" + str(new_anchor.snarl_id)) if extend_left else (str(new_anchor.snarl_id) + "-" + str(other_anchor.snarl_id))
                    new_anchor.add_snarl_id(new_snarl_id)
                    new_anchor.bp_occupied_start_node = (other_anchor.bp_occupied_start_node if extend_left else anchor.bp_occupied_start_node)
                    new_anchor.bp_occupied_end_node = (anchor.bp_occupied_end_node if extend_left else other_anchor.bp_occupied_end_node)
                    new_anchor.compute_bp_length()
                    new_anchor.reference_paths_covered = set(new_anchor.reference_paths_covered).intersection(set(other_anchor.reference_paths_covered))
                    # common_read_ids = set([read[READ_POSITION] for read in anchor.bp_matched_reads]).intersection(set([read[READ_POSITION] for read in other_anchor.bp_matched_reads]))
                    # find start, end of common reads
                    common_bp_matched_reads = {}
                    for read in new_anchor.bp_matched_reads:
                        if read[settings.READ_ID] in common_reads_ids:
                            common_bp_matched_reads[read[settings.READ_ID]] = read
                    for read in other_anchor.bp_matched_reads:
                        if read[settings.READ_ID] in common_reads_ids:
                            if settings.DEBUG:
                                print(f"processing read {read[settings.READ_ID]}", flush=True, file=stderr)
                            unpacked_read_id, unpacked_strand, unpacked_start, unpacked_end, unpacked_match_limit, unpacked_cs_left, unpacked_cs_right = common_bp_matched_reads[read[settings.READ_ID]]
                            if settings.DEBUG:
                                print(f"current anchor boundary before merge: {anchor!r} : {unpacked_start} - {unpacked_end}", flush=True, file=stderr)
                                print(f"other anchor boundary before merge: {other_anchor!r} : {read[settings.ANCHOR_START]} - {read[settings.ANCHOR_END]}", flush=True, file=stderr)
                            unpacked_start = min(read[settings.ANCHOR_START], unpacked_start)
                            unpacked_end = max(read[settings.ANCHOR_END], unpacked_end)
                            unpacked_cs_left = min(read[settings.CS_LEFT_AVAIL], unpacked_cs_left)
                            unpacked_cs_right = min(read[settings.CS_RIGHT_AVAIL], unpacked_cs_right)
                            common_bp_matched_reads[read[settings.READ_ID]] = [unpacked_read_id, unpacked_strand, unpacked_start, unpacked_end, unpacked_match_limit, unpacked_cs_left, unpacked_cs_right]
                            if settings.DEBUG:
                                print(f"DEBUG: In current snarl {current_snarl_id} after merging, read {read[settings.READ_ID]} has read info: {common_bp_matched_reads[read[settings.READ_ID]]}", flush=True, file=stderr)
                            if settings.DEBUG:
                                print(f"anchor boundary AFTER merge: {new_anchor!r} : {unpacked_start} - {unpacked_end}", flush=True, file=stderr)
                    
                    common_bp_matched_reads_list = list(common_bp_matched_reads.values())
                    new_anchor.bp_matched_reads = sorted(common_bp_matched_reads_list, key=lambda read: read[settings.READ_ID])    # set(anchor.bp_matched_reads).intersection(set(other_anchor.bp_matched_reads))
                    new_anchors_after_merging.append(new_anchor)

        # remove all anchors from both current and other snarls, as they are now replaced by new anchors having new snarl name
        for anchor in current_snarl_anchors:
            anchors_to_discard.add(f"{anchor!r}")
        for anchor in self.snarl_to_anchors_dictionary[other_snarl_id]:
            anchors_to_discard.add(f"{anchor!r}")
        # add the new snarl and its anchors in the snarl_to_anchors_dictionary
        new_snarl_id_after_merge = new_anchors_after_merging[0].snarl_id
        if settings.DEBUG:
            print(f"Inside snarl merging, have >=2 new anchors after merging {current_snarl_id} and {other_snarl_id}. New snarl id should be: {new_snarl_id_after_merge}", flush=True, file=stderr)
        self.snarl_to_anchors_dictionary[new_snarl_id_after_merge] = new_anchors_after_merging
        # updating the snarl_ids_sorted_list_up_to_date, and fixing it's iterator
        # inserting the new snarl id at the index iterator position
        snarl_ids_sorted_list_up_to_date.insert(snarl_ids_sorted_list_iterator_idx, new_snarl_id_after_merge)
        index_of_current_snarl_in_list = snarl_ids_sorted_list_up_to_date.index(current_snarl_id)
        index_of_other_snarl_in_list = snarl_ids_sorted_list_up_to_date.index(other_snarl_id)
        if index_of_current_snarl_in_list > index_of_other_snarl_in_list:
            snarl_ids_sorted_list_iterator_idx -= 1
        snarl_ids_sorted_list_up_to_date.remove(current_snarl_id)
        snarl_ids_sorted_list_up_to_date.remove(other_snarl_id)
        return (new_anchors_after_merging, snarl_ids_sorted_list_iterator_idx)
    

    def _extending_anchors_to_1_degree_node(self, node_handle_to_extend_to, current_snarl_boundary_handle, current_snarl_id, current_snarl_anchors, extend_left, anchors_to_discard, snarl_ids_sorted, snarl_ids_list_idx):
        """
        Extends anchors in a snarl into an adjacent 1-degree node. This function is called when we want to
        extend anchors into neighboring nodes that have only one connection (1-degree nodes).

        The function:
        1. Checks if there are enough reads that can be extended into the 1-degree node
        2. Extends the anchor boundaries to include the new node
        3. Updates read positions and coverage information
        4. Maintains anchor orientation and path information

        Parameters
        ----------
        node_handle_to_extend_to : node_handle
            Handle of the 1-degree node to extend into
        current_snarl_boundary_handle : node_handle
            Handle of the current snarl boundary node
        current_snarl_id : str
            ID of the current snarl being processed
        current_snarl_anchors : list
            List of Anchor objects in the current snarl
        extend_left : bool
            True if extending towards left, False if extending towards right
        anchors_to_discard : list
            List to store anchors that should be removed after extension
        snarl_ids_sorted : list
            List of all snarl IDs in sorted order
        snarl_ids_list_idx : int
            Current position in the snarl IDs list

        Returns
        -------
        list
            List of extended anchors if extension was successful, otherwise returns the original anchors
        """
        node_handle = node_handle_to_extend_to

        if settings.DEBUG:
            print(f"    ...currently checking extension in node ID {self.graph.get_id(node_handle)} which was extend_left = {extend_left} direction.")

        num_anchors_remaining_after_extension = 0
        anchors_to_extend = []
        for anchor in current_snarl_anchors:
            total_bp_matched_reads = len(anchor.bp_matched_reads)
            if settings.DEBUG:
                print(f"    ...anchor pre-extension is: {anchor!r}, total_bp_matched_reads = {total_bp_matched_reads}", flush=True, file=stderr)
            common_bp_matched_reads = []
            num_reads_that_can_be_extended = 0
            
            bp_added_upon_extention = (self.graph.get_length(current_snarl_boundary_handle) + 1) // 2 + (self.graph.get_length(node_handle_to_extend_to)) // 2 + 1
            # print(f"  bp_added_upon_extention = {bp_added_upon_extention}")

            for read in anchor.bp_matched_reads:
                read_id, read_strand, anchor_start, anchor_end, match_limit, cs_avail_left, cs_avail_right = read
                if settings.DEBUG:
                    print(f"    ...CHECKING for read {read_id}, read_strand = {read_strand}, anchor_start = {anchor_start}, anchor_end = {anchor_end}, match_limit = {match_limit}, cs_avail_left = {cs_avail_left}, cs_avail_right = {cs_avail_right}")
                # if read_strand == 0:    # FORWARD ORIENTATION
                    # check available bps in left direction
                if (
                    (read_strand == 0) and ((cs_avail_left if extend_left else cs_avail_right) >= bp_added_upon_extention)
                    or ((read_strand == 1) and ((cs_avail_right if extend_left else cs_avail_left) >= bp_added_upon_extention))
                ):
                    # print(f"    ...read 0 can be extended")
                    num_reads_that_can_be_extended += 1
                    # new anchor start position in read
                    read_start = (anchor_start - bp_added_upon_extention + 1) if extend_left else anchor_start
                    read_end =  anchor_end if extend_left else (anchor_end + bp_added_upon_extention - 1)
                    if read_strand == 0:
                        cs_in_left = (cs_avail_left - bp_added_upon_extention + 1) if extend_left else cs_avail_left
                        cs_in_right = cs_avail_right if extend_left else (cs_avail_right - bp_added_upon_extention + 1)
                    else:
                        cs_in_left = cs_avail_left if extend_left else (cs_avail_left - bp_added_upon_extention + 1)
                        cs_in_right = (cs_avail_right - bp_added_upon_extention + 1) if extend_left else cs_avail_right
                    # print(f"    ...APPENDING new anchor pos in read 0: read_start = {read_start}, read_end = {read_end}, cs_in_left = {cs_in_left}, cs_in_right = {cs_in_right}")
                    common_bp_matched_reads.append([read_id, read_strand, read_start, read_end, match_limit, cs_in_left, cs_in_right])
                                    
            # if more than threshold reads are dropped:
            if (num_reads_that_can_be_extended/total_bp_matched_reads >= settings.FRACTION_READS_FOR_SNARL_BOUNDARY_EXTENTION) and (num_reads_that_can_be_extended >= settings.MIN_READS_REQUIRED_FOR_BOUNDARY_EXTENSION):
                if settings.DEBUG:
                    print(f"    ...For anchor {anchor!r}, num_reads_that_can_be_extended = {num_reads_that_can_be_extended} and total_bp_matched_reads = {total_bp_matched_reads}, so anchor can be extended!")
                anchors_to_extend.append([anchor, common_bp_matched_reads])
                num_anchors_remaining_after_extension += 1

        # check if there are sufficient anchors remaining after extension
        if num_anchors_remaining_after_extension >= 2:
            # calculating the new snarl id
            new_snarl_id = current_snarl_id
            current_snarl_id = current_snarl_id.split("-")
            if extend_left:
                if "n" in current_snarl_id[0]:
                    node_number = int(current_snarl_id[0][1:]) + 1
                    new_snarl_id = "-".join(["n" + str(node_number)] + current_snarl_id[1:])
                else:
                    new_snarl_id = "-".join(["n1"] + current_snarl_id)
            else:
                if "n" in current_snarl_id[-1]:
                    node_number = int(current_snarl_id[-1][1:]) + 1
                    new_snarl_id = "-".join(current_snarl_id[:-1] + ["n" + str(node_number)])
                else:
                    new_snarl_id = "-".join(current_snarl_id + ["n1"])
            snarl_ids_sorted[snarl_ids_list_idx] = new_snarl_id

            anchors_list_to_extend = [anchor for anchor,_ in anchors_to_extend]
            for anchor in current_snarl_anchors:
                if anchor not in anchors_list_to_extend:
                    anchors_to_discard.append(anchor)
            
            # now do the extension
            for anchor, common_bp_matched_reads in anchors_to_extend:
                # do the extension here, in the same anchor object
                anchor.insert_node_through_extension(
                    Node(
                        self.graph.get_id(node_handle),
                        self.graph.get_length(node_handle),
                        not (self.graph.get_is_reverse(node_handle))
                    ), insert_left = extend_left
                )
                anchor.compute_bp_length()
                # if we extended on the left
                anchor.add_snarl_id(new_snarl_id)
                if settings.DEBUG:
                    print(f"    ....INSIDE: New snarl ID of anchor {anchor!r} after extending is {anchor.snarl_id}")
                anchor.bp_matched_reads = common_bp_matched_reads
            
            extended_anchors = [anchor for anchor, _ in anchors_to_extend]
            if settings.DEBUG:
                print(f"    ....INSIDE: extended anchors are {extended_anchors}")
            return extended_anchors
        
        else:
            if settings.DEBUG:
                print(" ....INSIDE: extension not done as num_anchors_remaining_after_extension < 2, so returned same old anchors = {current_snarl_anchors}")
            return current_snarl_anchors


    def _try_extension(self, current_snarl_anchors, current_snarl_id, other_snarl_id, anchors_to_discard, per_anchor_max_bps_to_extend, extend_left, extension_iteration):
        """
        Attempts to extend anchors in a snarl towards an adjacent snarl. This function handles the actual
        extension process, including:
        1. Calculating available base pairs for extension
        2. Updating anchor boundaries
        3. Adjusting read positions and coverage
        4. Handling read drops during extension

        Parameters
        ----------
        current_snarl_anchors : list
            List of Anchor objects in the current snarl
        current_snarl_id : str
            ID of the current snarl being processed
        other_snarl_id : str
            ID of the adjacent snarl to extend towards
        anchors_to_discard : list
            List to store anchors that should be removed after extension
        per_anchor_max_bps_to_extend : list
            List of maximum base pairs that can be extended for each anchor
        extend_left : bool
            True if extending towards left, False if extending towards right
        extension_iteration : int
            Current iteration number (0: no drops, 1: stricter drops, 2: relaxed drops)

        Returns
        -------
        list
            List of extended anchors if extension was successful, otherwise returns the original anchors
        """
        other_snarl_closest_node_id = 0
        an_other_snarl_anchor = self.snarl_to_anchors_dictionary[other_snarl_id][0]
        a_current_snarl_anchor = self.snarl_to_anchors_dictionary[current_snarl_id][0]

        if extend_left:
            other_snarl_closest_node_id =  max(an_other_snarl_anchor[0].id, an_other_snarl_anchor[-1].id)
            next_node_to_extend_node_id = current_snarl_boundary_node_id = min(a_current_snarl_anchor[0].id, a_current_snarl_anchor[-1].id)
        else:
            other_snarl_closest_node_id = min(an_other_snarl_anchor[0].id, an_other_snarl_anchor[-1].id)
            next_node_to_extend_node_id = current_snarl_boundary_node_id = max(a_current_snarl_anchor[0].id, a_current_snarl_anchor[-1].id)


        while True:
            cant_extend_more = False
            bp_occupied_current_snarl_boundary_node = a_current_snarl_anchor.bp_occupied_start_node if extend_left else a_current_snarl_anchor.bp_occupied_end_node
            bp_occupied_next_node = bp_occupied_current_snarl_boundary_node if (next_node_to_extend_node_id == current_snarl_boundary_node_id) else 0
            bp_occupied_other_snarl_boundary_node = an_other_snarl_anchor.bp_occupied_end_node if extend_left else an_other_snarl_anchor.bp_occupied_start_node
            next_node_handle = self.graph.get_handle(next_node_to_extend_node_id)
            if (
                next_node_to_extend_node_id == other_snarl_closest_node_id
            ):
                if settings.DEBUG:
                    print(f"bp_occupied_next_node = {bp_occupied_next_node}, bp_occupied_other_snarl_boundary_node = {bp_occupied_other_snarl_boundary_node}")
                bp_available_for_extension = self.graph.get_length(next_node_handle) - bp_occupied_next_node - bp_occupied_other_snarl_boundary_node     # bps available in current new node 
                if settings.DEBUG:
                    print(f"bp_available_for_extension = {bp_available_for_extension}")
                cant_extend_more = True
            else:    # means this node isn't a boundary node for another snarl, so at max, it can be consumed completely
                if settings.DEBUG:
                    print(f"bp_occupied_next_node = {bp_occupied_next_node}")
                bp_available_for_extension = self.graph.get_length(next_node_handle) - bp_occupied_next_node
                if settings.DEBUG:
                    print(f"bp_available_for_extension = {bp_available_for_extension}")

            # find overall bps to extend by doing min of bp_available_for_extension with all bps in per_anchor_max_bps_to_extend
            min_bp_among_cs_lines = min(per_anchor_max_bps_to_extend)
            if min_bp_among_cs_lines <= bp_available_for_extension:    # this means cs_avail will be completely consumed in this node 
                cant_extend_more = True
            final_bp_count_added_in_current_iteration = min(min_bp_among_cs_lines, bp_available_for_extension)

            # DOES a_current_snarl_anchor.basepairlength GET CORRECTLY UPDATED IN THIS WHILE LOOP? (yes)
            min_basepairlength_among_snarl_anchors = min(anchor.basepairlength for anchor in current_snarl_anchors)
            if (final_bp_count_added_in_current_iteration + min_basepairlength_among_snarl_anchors > settings.MIN_ANCHOR_LENGTH):
                final_bp_count_added_in_current_iteration = max(0, settings.MIN_ANCHOR_LENGTH - min_basepairlength_among_snarl_anchors)
                if settings.DEBUG:
                    print(f"We are in extension iteration {extension_iteration} and final_bp_count_added_in_current_iteration = {final_bp_count_added_in_current_iteration}")
                cant_extend_more = True

            if settings.DEBUG:
                print(f"    final_bp_count_added_in_current_iteration for snarl ID {current_snarl_id} = {final_bp_count_added_in_current_iteration}")
            
            # update per_anchor_max_bps_to_extend
            for idx in range(len(per_anchor_max_bps_to_extend)):
                per_anchor_max_bps_to_extend[idx] -= final_bp_count_added_in_current_iteration

            # do the extension
            if (next_node_to_extend_node_id != current_snarl_boundary_node_id) and (final_bp_count_added_in_current_iteration != 0):    # which means we need to first add this node to all anchors of this snarl
                # add this node to all anchors (in correct orientation)
                # get boundary node handle from id here
                for anchor in current_snarl_anchors:
                    anchor.insert_node_through_extension(
                        Node(
                            self.graph.get_id(next_node_handle),
                            self.graph.get_length(next_node_handle),
                            not (self.graph.get_is_reverse(next_node_handle))
                        ), insert_left = extend_left
                    )
                    if extend_left:
                        anchor.bp_occupied_start_node = 0
                    else:
                        anchor.bp_occupied_end_node = 0


            for anchor in current_snarl_anchors:
                # add bp_available_for_extension
                if extend_left:
                    anchor.bp_occupied_start_node += final_bp_count_added_in_current_iteration
                    if settings.DEBUG:
                        print(f"bp_occupied_start_node of anchor {anchor!r} is {anchor.bp_occupied_start_node}")
                else:
                    anchor.bp_occupied_end_node += final_bp_count_added_in_current_iteration
                    if settings.DEBUG:
                        print(f"bp_occupied_end_node of anchor {anchor!r} is {anchor.bp_occupied_end_node}")

                anchor.compute_bp_length()

                # update bp_matched_reads (also update cs_avail_left/right accordingly)
                new_bp_matched_reads = []
                for read in anchor.bp_matched_reads:
                    if extend_left:
                        path_orientation = (anchor._nodes[0].id < anchor._nodes[-1].id)
                        new_cs_avail_idx = (settings.CS_LEFT_AVAIL if ((read[settings.READ_STRAND] == 0 and path_orientation) or (read[settings.READ_STRAND] == 1 and not path_orientation)) else settings.CS_RIGHT_AVAIL)
                        
                        new_anchor_boundary_idx = settings.ANCHOR_START if path_orientation else settings.ANCHOR_END
                        if read[settings.READ_STRAND] == 1:
                            new_anchor_boundary_idx = (settings.ANCHOR_END if ((read[settings.READ_STRAND] == 0 and path_orientation) or (read[settings.READ_STRAND] == 1 and not path_orientation)) else settings.ANCHOR_START)
                        read[new_anchor_boundary_idx] = read[new_anchor_boundary_idx] - (final_bp_count_added_in_current_iteration if new_anchor_boundary_idx == settings.ANCHOR_START else (-final_bp_count_added_in_current_iteration))
                    else:
                        path_orientation = (anchor._nodes[0].id < anchor._nodes[-1].id)
                        new_cs_avail_idx = (settings.CS_RIGHT_AVAIL if ((read[settings.READ_STRAND] == 0 and path_orientation) or (read[settings.READ_STRAND] == 1 and not path_orientation)) else settings.CS_LEFT_AVAIL)
                        new_anchor_boundary_idx = settings.ANCHOR_END if path_orientation else settings.ANCHOR_START
                        if read[settings.READ_STRAND] == 1:
                            new_anchor_boundary_idx = (settings.ANCHOR_START if ((read[settings.READ_STRAND] == 0 and path_orientation) or (read[settings.READ_STRAND] == 1 and not path_orientation)) else settings.ANCHOR_END)
                        read[new_anchor_boundary_idx] = read[new_anchor_boundary_idx] - (final_bp_count_added_in_current_iteration if new_anchor_boundary_idx == settings.ANCHOR_START else (-final_bp_count_added_in_current_iteration))
                    
                    
                    new_cs_avail = read[new_cs_avail_idx] - final_bp_count_added_in_current_iteration
                    if new_cs_avail >= 0:
                        read[new_cs_avail_idx] = new_cs_avail
                        new_bp_matched_reads.append(read)
                    else:
                        if settings.OUTPUT_LOGGING_FILES:
                            if current_snarl_id not in self.anchor_read_tracking_dict:
                                self.anchor_read_tracking_dict[current_snarl_id] = dict()
                            if f"{anchor!r}" not in self.anchor_read_tracking_dict[current_snarl_id]:
                                self.anchor_read_tracking_dict[current_snarl_id][f"{anchor!r}"] = {}
                            if extension_iteration not in self.anchor_read_tracking_dict[current_snarl_id][f"{anchor!r}"]:
                                self.anchor_read_tracking_dict[current_snarl_id][f"{anchor!r}"][extension_iteration] = []
                            # add this read to anchor_read_tracking_dict
                            self.anchor_read_tracking_dict[current_snarl_id][f"{anchor!r}"][extension_iteration].append(read[settings.READ_ID])
                
                anchor.bp_matched_reads = new_bp_matched_reads

                # update anchor.compute_bp_length() to have correct calculation for boundary nodes
                # if anchor.basepairlength >= settings.MIN_ANCHOR_LENGTH:
                #     cant_extend_more = True

            # update current_snarl_boundary_node_id to next node to extend to, for next extension
            current_snarl_boundary_node_id = next_node_to_extend_node_id
            current_node_handle = next_node_handle
            # if flag to break is set, then break out of while loop here
            if cant_extend_more:
                break
            else:
                # follow edge in graph to get the next node id for the next iteration cycle
                current_node_out_degree = self.graph.get_degree(current_node_handle, extend_left)
                if current_node_out_degree == 1:    # extend here, not merge
                    self.graph.follow_edges(current_node_handle, extend_left, self.next_handle_iteratee)
                    next_node_handle = self.next_handle_expand_boundary
                    next_node_to_extend_node_id = self.graph.get_id(next_node_handle)
                else:
                    break
                
        return current_snarl_anchors


    def _get_max_read_drop(self, current_anchor_readcov: int, extension_iteration: int):
        """
        This function computes the allowed read drop for each anchor.
        
        """

        if extension_iteration == 1:
            ##### ITERATION-1: STRICTER READ-DROPS
            # How much more is the current anchor's read cov as compared to MIN_ANCHOR_READCOV
            allowed_read_drops = int(settings.DROP_FRACTION * current_anchor_readcov)
            if current_anchor_readcov <= settings.MIN_ANCHOR_READCOV:
                return 0
            elif current_anchor_readcov - allowed_read_drops <= settings.MIN_ANCHOR_READCOV:
                return current_anchor_readcov - settings.MIN_ANCHOR_READCOV
            else:
                return allowed_read_drops

        else:
            ##### ITERATION-2: RELAXED READ-DROPS, BUT PUSHES FOR TOUCHING MIN_ANCHOR_LENGTH
            return max(0, current_anchor_readcov - settings.MIN_ANCHOR_READCOV)

            
            # diff_readcov = current_anchor_readcov - settings.MIN_ANCHOR_READCOV

            # if diff_readcov < 0:
            #     return 0
            # elif diff_readcov < settings.MAX_READ_DROPS_ALLOWED:
            #     return diff_readcov
            # else:
            #     return settings.MAX_READ_DROPS_ALLOWED


    def _get_max_cs_avail_in_anchor(self, current_anchor_cs_avail_list: list, read_drops_allowed: int) -> int:
        """
        This function returns the max base pairs available for extension based on allowed read drops computed for that anchor
        """

        current_anchor_cs_avail_list_sorted = sorted(current_anchor_cs_avail_list)
        
        return(current_anchor_cs_avail_list_sorted[read_drops_allowed])

        
    def _extending_snarl_boundaries(self, current_snarl_anchors, current_snarl_id, snarl_ids_sorted, snarl_ids_list_idx, anchors_to_discard, extension_iteration):
        """
        Extends the boundaries of a snarl by attempting to extend its anchors. This function is the main
        coordinator for anchor extension, handling:
            1. Initial extension without read drops
            2. Extension with stricter read drop thresholds
            3. Extension with relaxed read drop thresholds
            4. Updating snarl boundaries and anchor information

        Parameters
        ----------
        current_snarl_anchors : list
            List of Anchor objects in the current snarl
        current_snarl_id : str
            ID of the current snarl being processed
        snarl_ids_sorted : list
            List of all snarl IDs in sorted order
        snarl_ids_list_idx : int
            Current position in the snarl IDs list
        anchors_to_discard : list
            List to store anchors that should be removed after extension
        extension_iteration : int
            Current iteration number (0: no drops, 1: stricter drops, 2: relaxed drops)

        Returns
        -------
        None
            The function modifies the anchors in place and updates the snarl boundaries
        """

        current_snarl_anchor_readcov = []    # storing read coverage of each anchor
        if settings.DEBUG:
            print(f"...extending left")
        for anchor in current_snarl_anchors:
            if settings.DEBUG:
                print(f"...current anchor's ({anchor!r}) basepairlength is {anchor.basepairlength} and bp_matched_reads are {anchor.bp_matched_reads}")
            current_snarl_anchor_readcov.append(len(anchor.bp_matched_reads))
        if settings.DEBUG:
            print(f"...current_snarl_anchor_readcov is {current_snarl_anchor_readcov}")

        allowed_read_drop_counts = [0] * len(current_snarl_anchors)    # stores for each anchor, how many reads can be dropped; 0 in the no_drop extension
        if extension_iteration != 0:
            allowed_read_drop_counts = [self._get_max_read_drop(current_anchor_readcov, extension_iteration) for current_anchor_readcov in current_snarl_anchor_readcov]
        if settings.DEBUG:
            print(f"...allowed_read_drop_counts is {allowed_read_drop_counts}")
        allowed_read_drop_counts_iterator = iter(allowed_read_drop_counts)
        per_anchor_max_bps_to_extend_left = []    # storing for each anchor, max cs_avail for extension in the left_direction 
        for anchor in current_snarl_anchors:    # for left extension
            # FIXME: This reliance on anchor node list order (to determine path orientation) is not robust. 
            # Later in extending_snarls_by_merging, we are liberally calling anchor.flip_anchors(). So, when using this comparison, disable merging. And later, revsiit that!!!!            per_read_cs_avail_list = [(read[settings.CS_LEFT_AVAIL] if (read[settings.READ_STRAND] == 0) else read[settings.CS_RIGHT_AVAIL]) for read in anchor.bp_matched_reads]
            path_orientation = (anchor._nodes[0].id < anchor._nodes[-1].id)
            per_read_cs_avail_list = [(read[settings.CS_LEFT_AVAIL] if ((read[settings.READ_STRAND] == 0 and path_orientation) or (read[settings.READ_STRAND] == 1 and not path_orientation)) else read[settings.CS_RIGHT_AVAIL]) for read in anchor.bp_matched_reads]
            
            for pos,read in enumerate(anchor.bp_matched_reads):
                if settings.DEBUG:
                    print(f"DEBUG: Extending left, for anchor {anchor!r}, read {read[settings.READ_ID]} has cs_avail on physical left side = {per_read_cs_avail_list[pos]}")
            
            current_allowed_read_drop_counts = next(allowed_read_drop_counts_iterator)
            per_anchor_max_bps_to_extend_left.append(self._get_max_cs_avail_in_anchor(per_read_cs_avail_list, current_allowed_read_drop_counts))
        if settings.DEBUG:
            print(f"...per_anchor_max_bps_to_extend_left is {per_anchor_max_bps_to_extend_left}")
        if snarl_ids_list_idx > 0:
            self._try_extension(current_snarl_anchors, current_snarl_id, snarl_ids_sorted[snarl_ids_list_idx - 1], anchors_to_discard, per_anchor_max_bps_to_extend_left, extend_left=True, extension_iteration = extension_iteration)   # for no_drop left extension
        
        if settings.DEBUG:
            print(f"...done extending left")
        
        for anchor in current_snarl_anchors:
            for pos,read in enumerate(anchor.bp_matched_reads):
                if settings.DEBUG:
                    path_orientation = (anchor._nodes[0].id < anchor._nodes[-1].id)
                    print(f"DEBUG: After physical left extension, for anchor {anchor!r}, read {read[settings.READ_ID]} has cs_avail on physical left side = {anchor.bp_matched_reads[pos][settings.CS_LEFT_AVAIL] if ((read[settings.READ_STRAND] == 0 and path_orientation) or (read[settings.READ_STRAND] == 1 and not path_orientation)) else anchor.bp_matched_reads[pos][settings.CS_RIGHT_AVAIL]}")


        for anchor in current_snarl_anchors:
            if settings.DEBUG:
                print(f"...current anchor's ({anchor!r}) new basepairlength is {anchor.basepairlength}, and new bp_matched_reads are {anchor.bp_matched_reads}")

        if settings.DEBUG:
            print(f"...extending right now")
        current_snarl_anchor_readcov = []    # storing read coverage of each anchor
        for anchor in current_snarl_anchors:
            current_snarl_anchor_readcov.append(len(anchor.bp_matched_reads))
        if settings.DEBUG:
            print(f"...current_snarl_anchor_readcov is {current_snarl_anchor_readcov}")

        allowed_read_drop_counts = [0] * len(current_snarl_anchors)    # stores for each anchor, how many reads can be dropped; 0 in the no_drop extension
        if extension_iteration != 0:
            allowed_read_drop_counts = [self._get_max_read_drop(current_anchor_readcov, extension_iteration) for current_anchor_readcov in current_snarl_anchor_readcov]
        allowed_read_drop_counts_iterator = iter(allowed_read_drop_counts)
        per_anchor_max_bps_to_extend_right = []    # storing for each anchor, max cs_avail for extension in the right_direction 
        for anchor in current_snarl_anchors:    # for right extension
            # FIXME: This reliance on anchor node list order (to determine path orientation) is not robust. 
            # Later in extending_snarls_by_merging, we are liberally calling anchor.flip_anchors(). So, when using this comparison, disable merging. And later, revsiit that!!!!            per_read_cs_avail_list = [(read[settings.CS_LEFT_AVAIL] if (read[settings.READ_STRAND] == 0) else read[settings.CS_RIGHT_AVAIL]) for read in anchor.bp_matched_reads]
            path_orientation = (anchor._nodes[0].id < anchor._nodes[-1].id)
            per_read_cs_avail_list = [(read[settings.CS_RIGHT_AVAIL] if ((read[settings.READ_STRAND] == 0 and path_orientation) or (read[settings.READ_STRAND] == 1 and not path_orientation)) else read[settings.CS_LEFT_AVAIL]) for read in anchor.bp_matched_reads]
            
            for pos,read in enumerate(anchor.bp_matched_reads):
                if settings.DEBUG:
                    print(f"DEBUG: Extending right, for anchor {anchor!r}, read {read[settings.READ_ID]} has cs_avail on physical right side = {per_read_cs_avail_list[pos]}")

            current_allowed_read_drop_counts = next(allowed_read_drop_counts_iterator)
            per_anchor_max_bps_to_extend_right.append(self._get_max_cs_avail_in_anchor(per_read_cs_avail_list, current_allowed_read_drop_counts))
        
        if snarl_ids_list_idx < len(snarl_ids_sorted) - 1:
            self._try_extension(current_snarl_anchors, current_snarl_id, snarl_ids_sorted[snarl_ids_list_idx + 1], anchors_to_discard, per_anchor_max_bps_to_extend_right, extend_left=False, extension_iteration = extension_iteration)   # for no_drop left extension
        if settings.DEBUG:
            print(f"...done extending right")
        
        for anchor in current_snarl_anchors:
            for pos,read in enumerate(anchor.bp_matched_reads):
                if settings.DEBUG:
                    path_orientation = (anchor._nodes[0].id < anchor._nodes[-1].id)
                    print(f"DEBUG: After physical right extension, for anchor {anchor!r}, read {read[settings.READ_ID]} has cs_avail on physical right side = {anchor.bp_matched_reads[pos][settings.CS_RIGHT_AVAIL] if ((read[settings.READ_STRAND] == 0 and path_orientation) or (read[settings.READ_STRAND] == 1 and not path_orientation)) else anchor.bp_matched_reads[pos][settings.CS_LEFT_AVAIL]}")

        
        for anchor in current_snarl_anchors:
            if settings.DEBUG:
                print(f"...current anchor's ({anchor!r}) new basepairlength is {anchor.basepairlength}, and new bp_matched_reads are {anchor.bp_matched_reads}")
            

    def _helper_extension_loop(self, snarl_ids_sorted, anchors_to_remove, extension_iteration, is_het_round=True):
        """
        This function is a helper for extending snarl boundaries in a loop. It is used to
        perform boundary extension for all snarls in a sorted list, allowing for read drops
        based on the specified extension iteration.
        Parameters
        ----------
        snarl_ids_sorted : list 
            List of all snarl IDs in sorted order
        anchors_to_remove : list
            List to store anchors that should be removed after extension
        extension_iteration : int
            Current iteration number (0: no drops, 1: stricter drops, 2: relaxed drops)
        is_het_round : bool
            Flag indicating whether to perform extension for het snarls (True) or hom snarls (False) (default is True)
        Returns
        -------
        None
            The function modifies the anchors in place and updates the snarl boundaries
        """
        snarl_ids_list_idx = 0 
        while snarl_ids_list_idx < len(snarl_ids_sorted):
            current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]
            if settings.DEBUG:
                print(f"Processing snarl ID: {current_snarl_id}")
            current_snarl_anchors = self.snarl_to_anchors_dictionary[current_snarl_id]
            if settings.DEBUG:
                print(f"..Running _extending_snarl_boundaries of snarl {current_snarl_id} with {len(current_snarl_anchors)} anchors")
            min_basepairlength_among_snarl_anchors = min([anchor.basepairlength for anchor in current_snarl_anchors])
            if min_basepairlength_among_snarl_anchors >= settings.MIN_ANCHOR_LENGTH:
                if settings.DEBUG:
                    print(f"Skipping snarl ID {current_snarl_id} as it's already sufficiently long. Length: {min_basepairlength_among_snarl_anchors}")
                snarl_ids_list_idx += 1
                continue
            if is_het_round:
                # we need to check if snarl is het snarl, i.e., has > 1 anchor
                if len(current_snarl_anchors) > 1:
                    self._extending_snarl_boundaries(current_snarl_anchors, current_snarl_id, snarl_ids_sorted, snarl_ids_list_idx, anchors_to_remove, extension_iteration)
            else:
                if len(current_snarl_anchors) == 1:
                    self._extending_snarl_boundaries(current_snarl_anchors, current_snarl_id, snarl_ids_sorted, snarl_ids_list_idx, anchors_to_remove, extension_iteration)
            if settings.DEBUG:
                print(f"...done extending snarl {current_snarl_id}")
            for anchor in self.snarl_to_anchors_dictionary[current_snarl_id]:
                if settings.DEBUG:
                    print(f"...anchor's ({anchor!r}) new basepairlength is {anchor.basepairlength}, and new bp_matched_reads are {anchor.bp_matched_reads}")
            snarl_ids_list_idx += 1


    def _helper_find_relevant_boundary_node_details_for_current_snarl(self, current_snarl_anchors, which_boundary):
        # calculate snarl max boundaries
        if which_boundary == 'left':
            min_left_node = 100000000000000
            left_bp_occupied = -1
            for anchor in current_snarl_anchors:
                left_node_idx_in_anchor_nodes_list = 0 if anchor._nodes[0].id < anchor._nodes[-1].id else -1
                if min_left_node > anchor._nodes[left_node_idx_in_anchor_nodes_list].id:
                    min_left_node = anchor._nodes[left_node_idx_in_anchor_nodes_list].id
                    left_bp_occupied = anchor.bp_occupied_start_node    # TODO: VERIFY IF bp_occupied_start_node is the correct value to use here, independent of the anchor orientation
                elif min_left_node == anchor._nodes[left_node_idx_in_anchor_nodes_list].id:
                    left_bp_occupied = max(left_bp_occupied, anchor.bp_occupied_start_node)
            return (min_left_node, left_bp_occupied)

        elif which_boundary == 'right':
            max_right_node = -1
            right_bp_occupied = -1
            for anchor in current_snarl_anchors:
                right_node_idx_in_anchor_nodes_list = 0 if anchor._nodes[0].id > anchor._nodes[-1].id else -1
                if max_right_node < anchor._nodes[right_node_idx_in_anchor_nodes_list].id:
                    max_right_node = anchor._nodes[right_node_idx_in_anchor_nodes_list].id
                    right_bp_occupied = anchor.bp_occupied_end_node    # TODO: VERIFY IF bp_occupied_end_node is the correct value to use here, independent of the anchor orientation
                elif max_right_node == anchor._nodes[right_node_idx_in_anchor_nodes_list].id:
                    right_bp_occupied = max(right_bp_occupied, anchor.bp_occupied_end_node)
            return (max_right_node, right_bp_occupied)

        else:
            raise ValueError("Invalid boundary type specified. Use 'left' or 'right'.")


    def _helper_find_bps_available_for_extension(self, current_snarl_id: int, other_snarl_id: Union[str, int], extend_left: bool):
        current_snarl_anchors_from_before_extension = self.before_extension_snarl_to_anchors_dictionary[current_snarl_id]
        current_snarl_boundary_node_id, current_snarl_boundary_node_bps_occupied = self._helper_find_relevant_boundary_node_details_for_current_snarl(current_snarl_anchors_from_before_extension, 'left' if extend_left else 'right')
        for anchor in current_snarl_anchors_from_before_extension:
            anchor.compute_bp_length()
        initial_bps_count = min([anchor.basepairlength for anchor in current_snarl_anchors_from_before_extension])
        current_bps_count = initial_bps_count
        other_snarl_anchors = self.snarl_to_anchors_dictionary[other_snarl_id]
        other_snarl_boundary_node_id, other_snarl_boundary_node_bp_occupied = self._helper_find_relevant_boundary_node_details_for_current_snarl(other_snarl_anchors, 'right' if extend_left else 'left')
        if settings.DEBUG:
            print(f"DEBUG: _helper_find_bps_available_for_extension - current_snarl_id: {current_snarl_id}, other_snarl_id: {other_snarl_id}, extend_left: {extend_left}")
            print(f"DEBUG: current_snarl_boundary_node_id: {current_snarl_boundary_node_id}, current_snarl_boundary_node_bps_occupied: {current_snarl_boundary_node_bps_occupied}")
            print(f"DEBUG: other_snarl_boundary_node_id: {other_snarl_boundary_node_id}, other_snarl_boundary_node_bp_occupied: {other_snarl_boundary_node_bp_occupied}")

        current_node_handle = self.graph.get_handle(current_snarl_boundary_node_id)

        while True:
            if current_bps_count >= settings.MIN_ANCHOR_LENGTH:
                return current_bps_count - initial_bps_count
            
            # just a safety check, never expected to hit this
            if (extend_left and (current_snarl_boundary_node_id < other_snarl_boundary_node_id)) or ((not extend_left) and (current_snarl_boundary_node_id > other_snarl_boundary_node_id)):
                break

            if settings.DEBUG:
                print(f"DEBUG: current_node_handle size: {self.graph.get_length(current_node_handle)}")
            node_length = self.graph.get_length(current_node_handle)
            current_occupied = current_snarl_boundary_node_bps_occupied
            other_occupied = other_snarl_boundary_node_bp_occupied if (current_snarl_boundary_node_id == other_snarl_boundary_node_id) else 0
            if settings.DEBUG:
                print(f"DEBUG: type(current_snarl_boundary_node_id): {type(current_snarl_boundary_node_id)}, type(other_snarl_boundary_node_id): {type(other_snarl_boundary_node_id)}")
            available_for_extension_bps_in_current_node = node_length - current_occupied - other_occupied
            if settings.DEBUG:
                print(f"DEBUG: Calculation: {node_length} - {current_occupied} - {other_occupied} = {available_for_extension_bps_in_current_node}")

            if available_for_extension_bps_in_current_node + current_bps_count >= settings.MIN_ANCHOR_LENGTH:
                if settings.DEBUG:
                    print(f"DEBUG: Enough base pairs available for extension in current node {current_snarl_boundary_node_id}. Returning {settings.MIN_ANCHOR_LENGTH - initial_bps_count}")
                return settings.MIN_ANCHOR_LENGTH - initial_bps_count
            elif current_snarl_boundary_node_id == other_snarl_boundary_node_id:
                result = available_for_extension_bps_in_current_node + current_bps_count - initial_bps_count
                if settings.DEBUG:
                    print(f"DEBUG: Returning {result} (available: {available_for_extension_bps_in_current_node}, current: {current_bps_count}, initial: {initial_bps_count})")
                return result
            else:
                current_bps_count += available_for_extension_bps_in_current_node
            
            # Find the next node for extension.
            # 1.) follow edge in graph to get the next node id for the next iteration cycle
            # 2.) update current_snarl_boundary_node_id and current_snarl_boundary_node_bps_occupied
            current_node_out_degree = self.graph.get_degree(current_node_handle, extend_left)
            if current_node_out_degree == 1:    # extend here, not merge
                self.graph.follow_edges(current_node_handle, extend_left, self.next_handle_iteratee)
                next_node_handle = self.next_handle_expand_boundary
                current_node_handle = next_node_handle
                current_snarl_boundary_node_id = self.graph.get_id(next_node_handle)
                current_snarl_boundary_node_bps_occupied = 0    # reset this for the next iteration
            else:
                if settings.DEBUG:
                    print(f"Stopping extension at snarl boundary node {current_snarl_boundary_node_id} with out-degree {current_node_out_degree} (branching point encountered)")
                return current_bps_count - initial_bps_count
                raise ValueError(f"Expected out-degree 1 for snarl boundary node {current_snarl_boundary_node_id}. Got {current_node_out_degree}")

        return current_bps_count - initial_bps_count


    def extend_and_insert_node(self, current_extended_anchor, current_snarl_boundary_node_id, current_snarl_boundary_bps_occupied, additional_bps_to_cover, extend_left):
        if settings.DEBUG:
            print(f"DEBUG: extend_and_insert_node called - node_id: {current_snarl_boundary_node_id}, bps_occupied: {current_snarl_boundary_bps_occupied}, additional_bps_to_cover: {additional_bps_to_cover}, extend_left: {extend_left}")
        current_node_handle = self.graph.get_handle(current_snarl_boundary_node_id)
        bps_extended_till_now = 0
        while bps_extended_till_now < additional_bps_to_cover:
            bps_available_current_node = self.graph.get_length(current_node_handle) - current_snarl_boundary_bps_occupied
            bps_to_extend_in_current_node = min(bps_available_current_node, additional_bps_to_cover - bps_extended_till_now)
            if settings.DEBUG:
                print(f"DEBUG: Node {current_snarl_boundary_node_id} - bps_available_current_node: {bps_available_current_node}, bps_to_extend_in_current_node: {bps_to_extend_in_current_node}, bps_extended_till_now: {bps_extended_till_now}")
            if bps_to_extend_in_current_node < 0:
                raise ValueError(f"Negative base pairs available for extension in node {current_snarl_boundary_node_id}. Check snarl boundary conditions.")
            # insert current node into the anchor if not already present
            current_anchor_boundary_node_to_compare = min(current_extended_anchor._nodes[0].id, current_extended_anchor._nodes[-1].id) if extend_left else max(current_extended_anchor._nodes[0].id, current_extended_anchor._nodes[-1].id)
            if settings.DEBUG:
                print(f"DEBUG: Anchor boundary node: {current_anchor_boundary_node_to_compare}, current node: {current_snarl_boundary_node_id}")
            if current_anchor_boundary_node_to_compare != current_snarl_boundary_node_id:
                if settings.DEBUG:
                    print(f"DEBUG: Inserting new node {current_snarl_boundary_node_id} into anchor")
                current_extended_anchor.insert_node_through_extension(
                    Node(
                        self.graph.get_id(current_node_handle),
                        self.graph.get_length(current_node_handle),
                        not (self.graph.get_is_reverse(current_node_handle))
                    ), insert_left=extend_left
                )
                # Update the occupied base pairs in the anchor
                if extend_left:
                    current_extended_anchor.bp_occupied_start_node = bps_to_extend_in_current_node
                else:
                    current_extended_anchor.bp_occupied_end_node = bps_to_extend_in_current_node
            else:
                if settings.DEBUG:
                    print(f"DEBUG: Updating existing node {current_snarl_boundary_node_id} in anchor")
                # If the node is already present in the anchor, we just need to update the occupied base pairs
                if extend_left:
                    current_extended_anchor.bp_occupied_start_node += bps_to_extend_in_current_node
                else:
                    current_extended_anchor.bp_occupied_end_node += bps_to_extend_in_current_node
            # Update the anchor's base pair length
            current_extended_anchor.compute_bp_length()
            if settings.DEBUG:
                print(f"DEBUG: Anchor length after update: {current_extended_anchor.basepairlength}")
            bps_extended_till_now += bps_to_extend_in_current_node
            if bps_extended_till_now >= additional_bps_to_cover:
                if settings.DEBUG:
                    print(f"DEBUG: Extension complete - reached target of {additional_bps_to_cover} bps")
                break
            # compute new node handle and node id for next iteration
            current_node_out_degree = self.graph.get_degree(current_node_handle, extend_left)
            if settings.DEBUG:
                print(f"DEBUG: Node {current_snarl_boundary_node_id} out-degree: {current_node_out_degree}")
            if current_node_out_degree == 1:
                if settings.DEBUG:
                    print(f"DEBUG: Following edge from node {current_snarl_boundary_node_id}")
                self.graph.follow_edges(current_node_handle, extend_left, self.next_handle_iteratee)
                next_node_handle = self.next_handle_expand_boundary
                current_node_handle = next_node_handle
                current_snarl_boundary_node_id = self.graph.get_id(current_node_handle)
                if settings.DEBUG:
                    print(f"DEBUG: Moved to next node: {current_snarl_boundary_node_id}")
                # Update the occupied base pairs in the anchor
                current_snarl_boundary_bps_occupied = 0    # reset this for the next iteration
            else:
                # Stop extension when we encounter a branching point (out-degree > 1)
                if settings.DEBUG:
                    print(f"DEBUG: Stopping extension at node {current_snarl_boundary_node_id} with out-degree {current_node_out_degree} (branching point encountered)")
                break

        if settings.DEBUG:
            print(f"DEBUG: extend_and_insert_node finished - total bps extended: {bps_extended_till_now}")
        return


    def update_current_anchor_details_with_new_boundary(self, current_extended_anchor: Anchor, current_nonextended_anchor: Anchor, best_subsequence_left_side_offset: int, best_subsequence_supporting_reads: list):
        current_extended_anchor.copy_from_anchor(current_nonextended_anchor)
        # Update the anchor's boundaries based on the best subsequence found
        # find bps to extend in the right direction here, before extending to the left direction
        bps_to_extend_left = best_subsequence_left_side_offset
        bps_to_extend_right = settings.MIN_ANCHOR_LENGTH - bps_to_extend_left - current_extended_anchor.basepairlength

        if bps_to_extend_left > 0:
            # We need to find the left boundary node details for extension
            current_snarl_boundary_node_id, current_snarl_boundary_bps_occupied = self._helper_find_relevant_boundary_node_details_for_current_snarl([current_extended_anchor], 'left')
            # Now extend to the right
            self.extend_and_insert_node(current_extended_anchor, current_snarl_boundary_node_id, current_snarl_boundary_bps_occupied, additional_bps_to_cover=bps_to_extend_left, extend_left=True)
        if bps_to_extend_right > 0:
            # We need to find the right boundary node details for extension
            current_snarl_boundary_node_id, current_snarl_boundary_bps_occupied = self._helper_find_relevant_boundary_node_details_for_current_snarl([current_extended_anchor], 'right')
            # Now extend to the right
            self.extend_and_insert_node(current_extended_anchor, current_snarl_boundary_node_id, current_snarl_boundary_bps_occupied, additional_bps_to_cover=bps_to_extend_right, extend_left=False)
        current_extended_anchor.bp_matched_reads = best_subsequence_supporting_reads
        return


    def extend_anchors_independently(self, snarl_ids_sorted: list, valid_anchors: list) -> list:
        """
        This function extends the anchors independently for each snarl.
        It first finds the bps available for extension in both directions, and then finds the best subsequence for each anchor.
        It then updates the anchor with the new boundaries and reads.
        """
        extension_round = ["HET", "HOM"]
        for round in extension_round:
            if settings.DEBUG:
                print(f"DEBUG: Processing {round}s for independent extension")
            for current_snarl_idx, current_snarl_id in enumerate(snarl_ids_sorted):
                if (len(self.snarl_to_anchors_dictionary[current_snarl_id]) == 1 and round == "HOM") or (len(self.snarl_to_anchors_dictionary[current_snarl_id]) > 1 and round == "HET"):
                    if settings.DEBUG:
                        print(f"DEBUG: Processing snarl {current_snarl_id} for independent extension (index {current_snarl_idx})")
                    if isinstance(current_snarl_id, str) and ('-' in current_snarl_id):  # merged snarl
                        if settings.DEBUG:
                            print(f"DEBUG: Skipping merged snarl {current_snarl_id} for independent extension")
                        continue
                    if len(self.snarl_to_anchors_dictionary[current_snarl_id]) == 1:
                        if settings.DEBUG:
                            print(f"DEBUG: Skipping homozygous snarl {current_snarl_id} for independent extension")
                        continue

                    should_consider_for_extension = False
                    for anchor in self.snarl_to_anchors_dictionary[current_snarl_id]:
                        anchor.compute_bp_length()
                        if anchor.basepairlength < settings.MIN_ANCHOR_LENGTH:
                            should_consider_for_extension = True
                    if not should_consider_for_extension:
                        if settings.DEBUG:
                            print(f"DEBUG: Skipping snarl {current_snarl_id} for independent extension - all anchors are already long enough")
                        continue
                    
                    # find bps available to extend in both directions
                    left_snarl_idx = (current_snarl_idx - 1) if (current_snarl_idx > 0) else -1
                    if left_snarl_idx < 0:
                        bps_available_for_extension_on_left_side = 0
                    else:
                        bps_available_for_extension_on_left_side = self._helper_find_bps_available_for_extension(current_snarl_id, snarl_ids_sorted[left_snarl_idx], extend_left=True) if (left_snarl_idx >= 0) else 0
                        if settings.DEBUG:
                            print(f"DEBUG: Snarl {current_snarl_id} - left_snarl: {snarl_ids_sorted[left_snarl_idx]}, bps_available_for_extension_on_left_side: {bps_available_for_extension_on_left_side}")
                    
                    right_snarl_idx = (current_snarl_idx + 1) if (current_snarl_idx < len(snarl_ids_sorted) - 1) else len(snarl_ids_sorted)
                    if right_snarl_idx >= len(snarl_ids_sorted):
                        bps_available_for_extension_on_right_side = 0
                    else:
                        bps_available_for_extension_on_right_side = self._helper_find_bps_available_for_extension(current_snarl_id, snarl_ids_sorted[right_snarl_idx], extend_left=False) if (right_snarl_idx < len(snarl_ids_sorted)) else 0
                        if settings.DEBUG:
                            print(f"DEBUG: Snarl {current_snarl_id} - right_snarl: {snarl_ids_sorted[right_snarl_idx]}, bps_available_for_extension_on_right_side: {bps_available_for_extension_on_right_side}")
                    
                    min_anchor_length = min([anchor.basepairlength for anchor in self.before_extension_snarl_to_anchors_dictionary[current_snarl_id]])
                    total_available = bps_available_for_extension_on_left_side + bps_available_for_extension_on_right_side + min_anchor_length
                    if settings.DEBUG:
                        print(f"DEBUG: Snarl {current_snarl_id} - min_anchor_length: {min_anchor_length}, total_available: {total_available}, MIN_ANCHOR_LENGTH: {settings.MIN_ANCHOR_LENGTH}")
                    
                    if total_available < settings.MIN_ANCHOR_LENGTH:
                        if settings.DEBUG:
                            print(f"DEBUG: Skipping snarl {current_snarl_id} - insufficient total available bps")
                        continue

                    if settings.DEBUG:
                        print(f"DEBUG: Processing snarl {current_snarl_id} for independent extension")
                    for current_anchor_idx, current_anchor in enumerate(self.before_extension_snarl_to_anchors_dictionary[current_snarl_id]):
                        if settings.DEBUG:
                            print(f"DEBUG: Processing anchor {current_anchor!r} for snarl {current_snarl_id}")
                        # record the offset of the best subsequence (i.e., the one with most reads retained) STARTING FROM THE current snarl left boundary start node (the one before any kind of extension), and increasing in the LEFT DIRECTION            
                        best_subsequence_left_side_offset = settings.MIN_ANCHOR_LENGTH
                        best_subsequence_supporting_reads = []
                        current_anchor.compute_bp_length()
                        for current_subsequence_left_side_offset in range(max(0, min(settings.MIN_ANCHOR_LENGTH - current_anchor.basepairlength, bps_available_for_extension_on_left_side)), max(0, settings.MIN_ANCHOR_LENGTH - (bps_available_for_extension_on_right_side + current_anchor.basepairlength)), -1):
                            reads_supporting_current_subsequence = []
                            for read in current_anchor.bp_matched_reads:
                                right_side_cs_avail_required = settings.MIN_ANCHOR_LENGTH - current_subsequence_left_side_offset - current_anchor.basepairlength
                                path_orientation = (current_anchor[0].id < current_anchor[-1].id)
                                if (path_orientation and (read[settings.READ_STRAND] == 0)) or ((not path_orientation) and (read[settings.READ_STRAND] == 1)):  # forward strand
                                    if read[settings.CS_LEFT_AVAIL] >= current_subsequence_left_side_offset and read[settings.CS_RIGHT_AVAIL] >= right_side_cs_avail_required:
                                        reads_supporting_current_subsequence.append(read)
                                else:  # reverse strand
                                    # reverse strand reads have their left and right CS avail swapped
                                    # so we check CS_RIGHT_AVAIL for left side offset and CS_LEFT_AVAIL for right side offset
                                    if read[settings.CS_RIGHT_AVAIL] >= current_subsequence_left_side_offset and read[settings.CS_LEFT_AVAIL] >= right_side_cs_avail_required:
                                        reads_supporting_current_subsequence.append(read)
                                
                            if len(reads_supporting_current_subsequence) > len(best_subsequence_supporting_reads):
                                best_subsequence_left_side_offset = current_subsequence_left_side_offset
                                best_subsequence_supporting_reads = reads_supporting_current_subsequence
                        
                        # now update the current_anchor to have the boundaries defined by the best_subsequence_left_side_offset, and reads as best_subsequence_supporting_reads
                        if len(best_subsequence_supporting_reads) < settings.MIN_ANCHOR_READCOV_FOR_INDEPENDENT_ANCHOR_EXTENSION:
                            if settings.DEBUG:
                                print(f"DEBUG: Skipping anchor {current_anchor!r} for snarl {current_snarl_id} - insufficient read coverage for independent extension ({len(best_subsequence_supporting_reads)} < {settings.MIN_ANCHOR_READCOV_FOR_INDEPENDENT_ANCHOR_EXTENSION})")
                            # this means we cannot extend this anchor as it will be too short
                            continue

                        if settings.OUTPUT_LOGGING_FILES:
                            if current_snarl_id not in self.independent_anchor_extension_tracking_dict:
                                self.independent_anchor_extension_tracking_dict[current_snarl_id] = dict()
                            if current_anchor_idx not in self.independent_anchor_extension_tracking_dict[current_snarl_id]:
                                self.independent_anchor_extension_tracking_dict[current_snarl_id][current_anchor_idx] = dict()
                            self.independent_anchor_extension_tracking_dict[current_snarl_id][current_anchor_idx]["primary_anchor"] = [f"{self.before_extension_snarl_to_anchors_dictionary[current_snarl_id][current_anchor_idx]!r}", {"anchor_length": self.before_extension_snarl_to_anchors_dictionary[current_snarl_id][current_anchor_idx].basepairlength}, {"read_cov": len(self.before_extension_snarl_to_anchors_dictionary[current_snarl_id][current_anchor_idx].bp_matched_reads)}]
                            self.independent_anchor_extension_tracking_dict[current_snarl_id][current_anchor_idx]["extension_around_sentinel"] = [f"{self.snarl_to_anchors_dictionary[current_snarl_id][current_anchor_idx]!r}", {"anchor_length": self.snarl_to_anchors_dictionary[current_snarl_id][current_anchor_idx].basepairlength}, {"read_cov": len(self.snarl_to_anchors_dictionary[current_snarl_id][current_anchor_idx].bp_matched_reads)}]

                        # calculate correct boundaries of the current anchor in the reads belonging to best_subsequence_supporting_reads, and also their cs_avails
                        for read_idx, read in enumerate(best_subsequence_supporting_reads):
                            path_orientation = (current_anchor[0].id < current_anchor[-1].id)
                            ### calculating left and right anchor boundary indices in the read, based on the new definition of read strand 
                            # previously (wrong):
                            # read[settings.ANCHOR_START] -= best_subsequence_left_side_offset
                            new_left_anchor_boundary_idx = settings.ANCHOR_START if path_orientation else settings.ANCHOR_END
                            if read[settings.READ_STRAND] == 1:
                                new_left_anchor_boundary_idx = (settings.ANCHOR_END if (read[settings.READ_STRAND] == 1 and not path_orientation) else settings.ANCHOR_START)
                            read[new_left_anchor_boundary_idx] = read[new_left_anchor_boundary_idx] - (best_subsequence_left_side_offset if new_left_anchor_boundary_idx == settings.ANCHOR_START else (-best_subsequence_left_side_offset))

                            # previously (wrong):
                            # read[settings.ANCHOR_END] += settings.MIN_ANCHOR_LENGTH - best_subsequence_left_side_offset - current_anchor.basepairlength
                            best_subsequence_right_side_offset = settings.MIN_ANCHOR_LENGTH - best_subsequence_left_side_offset - current_anchor.basepairlength
                            new_right_anchor_boundary_idx = settings.ANCHOR_END if path_orientation else settings.ANCHOR_START
                            if read[settings.READ_STRAND] == 1:
                                new_right_anchor_boundary_idx = (settings.ANCHOR_START if (read[settings.READ_STRAND] == 1 and not path_orientation) else settings.ANCHOR_END)
                            read[new_right_anchor_boundary_idx] = read[new_right_anchor_boundary_idx] - (best_subsequence_right_side_offset if new_right_anchor_boundary_idx == settings.ANCHOR_START else (-best_subsequence_right_side_offset))
                            
                            # cs left and right avail calculations
                            if (path_orientation and (read[settings.READ_STRAND] == 0)) or ((not path_orientation) and (read[settings.READ_STRAND] == 1)):
                                read[settings.CS_LEFT_AVAIL] -= best_subsequence_left_side_offset
                                read[settings.CS_RIGHT_AVAIL] -= (settings.MIN_ANCHOR_LENGTH - best_subsequence_left_side_offset - current_anchor.basepairlength)
                            else:
                                read[settings.CS_RIGHT_AVAIL] -= best_subsequence_left_side_offset
                                read[settings.CS_LEFT_AVAIL] -= (settings.MIN_ANCHOR_LENGTH - best_subsequence_left_side_offset - current_anchor.basepairlength)
                            best_subsequence_supporting_reads[read_idx] = read

                        # updates the original anchor (i.e., the instance which had been extended previously through drops) with the new boundaries
                        # this way, we don't have to create a new anchor object and worry about managing its presence in valid_anchors.
                        if settings.DEBUG:
                            print(f"Selected best subsequence left side offset: {best_subsequence_left_side_offset} for snarl {current_snarl_id} anchor {current_anchor!r}")
                        self.update_current_anchor_details_with_new_boundary(self.snarl_to_anchors_dictionary[current_snarl_id][current_anchor_idx], self.before_extension_snarl_to_anchors_dictionary[current_snarl_id][current_anchor_idx], best_subsequence_left_side_offset, best_subsequence_supporting_reads)
                        if settings.OUTPUT_LOGGING_FILES:
                            self.independent_anchor_extension_tracking_dict[current_snarl_id][current_anchor_idx]["fake_anchor_generation"] = [f"{self.snarl_to_anchors_dictionary[current_snarl_id][current_anchor_idx]!r}", {"anchor_length": self.snarl_to_anchors_dictionary[current_snarl_id][current_anchor_idx].basepairlength}, {"read_cov": len(self.snarl_to_anchors_dictionary[current_snarl_id][current_anchor_idx].bp_matched_reads)}]

        return valid_anchors
    

    def _helper_determine_if_snarl_underwent_independent_extension(self, current_snarl_anchors: list) -> bool:

        #check if current snarl underwent independent extension
        last_anchor_start = min(current_snarl_anchors[0][0].id, current_snarl_anchors[0][-1].id)
        last_anchor_end = max(current_snarl_anchors[0][0].id, current_snarl_anchors[0][-1].id)
        last_anchor_bp_occupied_start, last_anchor_bp_occupied_end = current_snarl_anchors[0].bp_occupied_start_node, current_snarl_anchors[0].bp_occupied_end_node
        for anchor in current_snarl_anchors[1:]:
            if [min(anchor[0].id, anchor[-1].id), max(anchor[0].id, anchor[-1].id), anchor.bp_occupied_start_node, anchor.bp_occupied_end_node] == [last_anchor_start, last_anchor_end, last_anchor_bp_occupied_start, last_anchor_bp_occupied_end]:
                continue
            else:
                return True
        return False
    

    def _helper_fetch_list_of_anchor_sequences(self, current_snarl_anchors: list) -> list:
        """
        This function extracts the sequence of the anchor from the fasta file.
        It also reverses the sequence if the read is on the reverse strand.
        It returns a list of sequences, one for each anchor.
        """
        current_snarl_anchors_sequence_list = []
        # translating all reads into READ_STRAND=0
        for anchor in current_snarl_anchors:
            read = anchor.bp_matched_reads[0] # Extracting the first read from the anchor to get the sequence
            read_seq = helpers.extract_sequence(fasta_file=self.fasta_path, read_id=read[settings.READ_ID])
            if read_seq is not None:
                # Can the extracted read_seq from the fasta file and this read entry have opposite strands?
                anchor_slice_in_read = ""
                if read[settings.READ_STRAND] == 0:
                    anchor_slice_in_read = read_seq[read[settings.ANCHOR_START]:read[settings.ANCHOR_END]]
                else:
                    anchor_slice_in_read = read_seq[::-1][read[settings.ANCHOR_START]:read[settings.ANCHOR_END]]  # reverse the slice if the read is on the reverse strand
                    anchor_slice_in_read = helpers.complement(anchor_slice_in_read)
                current_snarl_anchors_sequence_list.append(anchor_slice_in_read)
            else:
                raise ValueError(f"Read {read[settings.READ_ID]} not found in fasta file")
        return current_snarl_anchors_sequence_list
    

    def _helper_extract_canonical_signature(self, anchor_seq: str, repeat_segments_offsets_list: list) -> tuple:
        """
        This function extracts the canonical signature of an anchor sequence, which is the sequence of the anchor sequence
        that is not part of the computed repeat segment(s) found by `sdust`.
        """
        repeat_segments_offsets_list.sort()
        canonical_signature_list = []
        last_end_offset = 0
        for repeat_segment in repeat_segments_offsets_list:
            if repeat_segment[0] > last_end_offset:
                canonical_signature_list.append(anchor_seq[last_end_offset:repeat_segment[0]])
            last_end_offset = repeat_segment[1]
        if len(anchor_seq) > last_end_offset:
            canonical_signature_list.append(anchor_seq[last_end_offset:])
        return tuple(canonical_signature_list)


    def _helper_find_low_complexity_regions(self, anchor_seq: str, w: int, t: int) -> list:
        """
        This function uses the `sdust` tool (https://github.com/lh3/sdust) to find low-complexity region ranges in a sequence

        Args:
            anchor_seq (str): The anchor sequence to find low-complexity regions in
            w (int): The window size for the sdust tool
            t (int): The threshold for the sdust tool

        Returns:
            list: A list of tuples, each containing the start and end indices of the low-complexity region
            in the anchor sequence
            Note: If no repeats are found, the function returns an empty list
        """
        
        # Save anchor sequence to a tmp fasta file
        temp_file = None
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fa") as temp_f:
            temp_file = temp_f.name
            temp_f.write(f">temp_seq\n{anchor_seq}\n")
        
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        sdust_executable_path = os.path.join(project_root, "bin", "sdust")
        command = [sdust_executable_path, "-w", str(w), "-t", str(t), temp_file]

        process = subprocess.run(command, capture_output=True, text=True, check=True)

        intervals = []
        for line in process.stdout.strip().split('\n'):  # each line is a repeat segment within the anchor sequence
            interval_parts = line.split()
            if len(interval_parts) != 0:
                repeat_tuple = (int(interval_parts[1]), int(interval_parts[2]))
                intervals.append(repeat_tuple)
        os.remove(temp_file)
        return intervals


    def prune_repeat_anchors(self, snarl_ids_sorted: list, valid_anchors: list) -> list:
        """
        This function prunes repeat anchors from the valid_anchors list.
        """
        if settings.DEBUG:
            print(f"#### PRUNING REPEAT ANCHORS ######")
        valid_anchors_after_pruning = []
        anchors_pruned = []
        snarl_id_idx = 0
        while snarl_id_idx < len(snarl_ids_sorted):
            snarl_id = snarl_ids_sorted[snarl_id_idx]
            if settings.DEBUG:
                print()
                print(f"#### PRUNING: Processing snarl {snarl_id} ######")
            # not considering merged anchors, directly adding them to final valid anchors list
            if isinstance(snarl_id, str) and '-' in snarl_id:
                if settings.DEBUG:
                    print(f"#### PRUNING: Snarl {snarl_id} is a merged snarl. Skipping... ######")
                [valid_anchors_after_pruning.append(anchor) for anchor in self.snarl_to_anchors_dictionary[snarl_id]]
                snarl_id_idx += 1
                continue
            current_snarl_anchors = self.snarl_to_anchors_dictionary[snarl_id]
            if len(current_snarl_anchors) == 1:  # homozygous snarl
                if settings.DEBUG:
                    print(f"#### PRUNING: Snarl {snarl_id} is a homozygous snarl. Skipping... ######")
                valid_anchors_after_pruning.append(current_snarl_anchors[0])
                snarl_id_idx += 1
                continue
            
            if self._helper_determine_if_snarl_underwent_independent_extension(current_snarl_anchors):
                if settings.DEBUG:
                    print(f"#### PRUNING: Snarl {snarl_id} was extended in independent extension. Skipping... ######")
                # means this snarl was extended in independent extension
                [valid_anchors_after_pruning.append(anchor) for anchor in current_snarl_anchors]
                snarl_id_idx += 1
                continue
            else:
                # determine if this snarl has repeat anchors
                found_repeat_anchor = False
                non_repeat_sequences_count_dict = {}
                current_snarl_anchors_sequence_list = self._helper_fetch_list_of_anchor_sequences(current_snarl_anchors)
                for anchor_idx, anchor_seq in enumerate(current_snarl_anchors_sequence_list):
                    if settings.DEBUG:
                        print(f"#### PRUNING: Processing anchor: {current_snarl_anchors[anchor_idx]!r}, sequence: {anchor_seq} ######")
                    repeat_segments_offsets_list = self._helper_find_low_complexity_regions(anchor_seq, w=40, t=4)
                    canonical_signature_tuple = self._helper_extract_canonical_signature(anchor_seq, repeat_segments_offsets_list)
                    if settings.DEBUG:
                        print(f"#### PRUNING: Canonical signature: {canonical_signature_tuple} ######")
                    non_repeat_sequences_count_dict[canonical_signature_tuple] = non_repeat_sequences_count_dict.get(canonical_signature_tuple, 0) + 1
                    if non_repeat_sequences_count_dict[canonical_signature_tuple] > 1:
                        found_repeat_anchor = True
                        break
                if not found_repeat_anchor:
                    if settings.DEBUG:
                        print(f"#### PRUNING: Snarl {snarl_id} does not have repeat anchors. Skipping... ######")
                    [valid_anchors_after_pruning.append(anchor) for anchor in current_snarl_anchors]
                else:
                    if settings.DEBUG:
                        print(f"#### PRUNING: Snarl {snarl_id} has repeat anchors. Adding all anchors to anchors_pruned, and removing snarl from snarl_ids_sorted... ######")
                    snarl_ids_sorted.remove(snarl_id)
                    snarl_id_idx -= 1
                    [anchors_pruned.append(anchor) for anchor in current_snarl_anchors]
                snarl_id_idx += 1
        return valid_anchors_after_pruning, anchors_pruned


    def extend_and_merge_snarls(self, valid_anchors: list) -> list:
        """
        This function performs snarl boundary extension and merging
        """
        anchors_to_remove = set()   # {(snarl_id, anchor)}

        # NOTE: This step is adding to the overhead the most.
        self.before_extension_snarl_to_anchors_dictionary = copy.deepcopy(self.snarl_to_anchors_dictionary)
        
        ### First, performing perfect bp match extension (no read drop allowed) for all snarls
        if settings.DEBUG:
            print(f"#### RUNNING EXTENSION WITH NO DROPS, FRACTIONAL ALLOWED DROPS AND THEN MORE DROPS FOR HET ANCHORS ONLY ####")
        t0 = time.time()
        self._helper_extension_loop(self.snarl_ids_sorted, anchors_to_remove, extension_iteration=0, is_het_round=True)
        self._helper_extension_loop(self.snarl_ids_sorted, anchors_to_remove, extension_iteration=1, is_het_round=True)
        self._helper_extension_loop(self.snarl_ids_sorted, anchors_to_remove, extension_iteration=2, is_het_round=True)
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(f"..Regular extension of HET anchors took {time.time() - t0} seconds", flush=True, file=stderr)

        if settings.DEBUG:
            print(f"#### RUNNING EXTENSION WITH NO DROPS, FRACTIONAL ALLOWED DROPS AND THEN MORE DROPS FOR HOM ANCHORS ONLY ####")
        t1 = time.time()
        self._helper_extension_loop(self.snarl_ids_sorted, anchors_to_remove, extension_iteration=0, is_het_round=False)
        self._helper_extension_loop(self.snarl_ids_sorted, anchors_to_remove, extension_iteration=1, is_het_round=False)
        self._helper_extension_loop(self.snarl_ids_sorted, anchors_to_remove, extension_iteration=2, is_het_round=False)
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(f"..Regular extension of HOM anchors took {time.time() - t1} seconds", flush=True, file=stderr)

        if settings.DEBUG:
            print(f"#### TRY TO MERGE SHORTER ANCHORS ######")
        t2 = time.time()
        valid_anchors = self.merge_anchors(valid_anchors, anchors_to_remove, self.snarl_ids_sorted, merging_round=0)
        
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(f"..Merging anchors took {time.time() - t2} seconds", flush=True, file=stderr)

        t3 = time.time()
        if settings.DEBUG:
            print(f"#### RUNNING INDEPENDENT ANCHOR EXTENSION ######")
        # Note: Now that snarl boundaries will not be the same as its anchors' boundaries, we will use 
        # self._helper_find_relevant_boundary_node_details_for_current_snarl() to calculate snarl's extreme boundaries on the fly
        
        valid_anchors = self.extend_anchors_independently(snarl_ids_sorted=self.snarl_ids_sorted, valid_anchors=valid_anchors)
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(f"..Independent anchor extension took {time.time() - t3} seconds", flush=True, file=stderr)


        # Note: valid_anchors is a list of lists, where each nested list contains an anchor object and a list of reads
        for idx in range(len(valid_anchors)):
            anchor = valid_anchors[idx][0]
            reads = [read[:4] for read in anchor.bp_matched_reads]
            valid_anchors[idx][1] = reads
                
        return valid_anchors    #### change this later to calculate valid_anchors_extended, when we will have anchor drops because of merging


    def merge_anchors(self, valid_anchors: list, anchors_to_remove: set, snarl_ids_sorted: list, merging_round: int) -> list:
        """
        * Iterate over shorter anchors, find adjacent snarls (+1/-1). If read drop from one snarl to the other is within the defined threshold,
        then merge the snarls. Get all combinations of anchors (required it belongs to atleast one path) and re-define this as a new anchor,
        with common reads and newly computed bplength, snarl_id and pathnames.

        * If both left, right anchors are available for merging, choose the direction where read coverage drop is minimum
            # get snarl direction (0 -> forward, 1 -> backward)
            # when merging anchors, if orientation of other anchor is opposite to that of current, do the following:-
            # 1.) flip other anchor
            # 2.) if other anchor can be found by extending in same direction as current anchor, then add other anchor to the end of current anchor. else, add to beginning.

        """
        valid_anchor_extended = []
        # newsnarls_to_anchors_dictionary = {}
        # anchors_to_remove = []   # {(snarl_id, anchor)}
        snarl_orientation = True

        # snarl_ids_sorted = sorted(list(self.snarl_to_anchors_dictionary.keys()))
        snarl_ids_list_idx = 0    # we need to use this index counter (and can't simply use an iterator) because we will be inserting/deleting the snarl_ids_sorted list on the go
        # Note: Remember to update the snarl_ids_list_idx appropriately when merging snarls

        # FIXME: Instead of iterating overall all snarls to find the ones to merge, we could use a shorter list of short snarls, which could be populated during extension.
        while snarl_ids_list_idx < len(snarl_ids_sorted):
            # print(f"Current snarl id: {snarl_ids_sorted[snarl_ids_list_idx]}", flush=True, file=stderr)
            current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]
            # current_snarl_first_snarl_id = current_snarl_id if isinstance(current_snarl_id, int) else int(current_snarl_id.split('-')[0])
            # if current_snarl_first_snarl_id < 39045:
            #     snarl_ids_list_idx += 1
            #     continue
            # if current_snarl_first_snarl_id > 39392:
            #     break
            last_snarl_id = -1
            current_snarl_anchors = self.snarl_to_anchors_dictionary[current_snarl_id]
            min_anchor_length_in_snarl = min([anchor.basepairlength for anchor in current_snarl_anchors])
            # if len(current_snarl_anchors) > 1 and min_anchor_length_in_snarl < settings.MIN_ANCHOR_LENGTH:

            while (
                len(current_snarl_anchors) > 1 
                and (min_anchor_length_in_snarl < settings.MIN_ANCHOR_LENGTH
                and last_snarl_id != current_snarl_id)
            ):
                current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]
                last_snarl_id = current_snarl_id

                # recalculating left and right nodes in graph for extension/merging
                if snarl_orientation:
                    # even if 0-th anchor of snarl is reversed, snarl_start and snarl_end should be in increasing order of node ids (if snarl_orientation = True)
                    if settings.DEBUG:
                        print(f"########################", end="\\n", flush=True, file=stderr)
                        print(f"Processing SNARL ID: {current_snarl_id}", flush=True, file=stderr)
                    current_snarl_start_id = min(current_snarl_anchors[0][0].id, current_snarl_anchors[0][-1].id)   # change this variable to current_snarl_start_node_id
                    if settings.DEBUG:
                        print(f"Current snarl start node is {current_snarl_start_id}", flush=True, file=stderr)
                    go_left_bool = True
                    current_snarl_start_handle = self.graph.get_handle(current_snarl_start_id)
                    extend_left = True

                    # Get left snarl's ID
                    # FIXME: We don't know what left of current snarl is. We need to infer this in some other way.
                    left_snarl_id = snarl_ids_sorted[snarl_ids_sorted.index(current_snarl_id) - 1] if (snarl_ids_sorted.index(current_snarl_id) - 1 >= 0) else -1 
                    left_snarl_end_node_id = -1
                    if (
                            self.snarl_to_anchors_dictionary.get(left_snarl_id) != None
                            and len(self.snarl_to_anchors_dictionary[left_snarl_id]) > 0
                        ):
                        # Get left snarl's end node. We need this to later check if our current snarl could be extended in the left direction 
                        # or not, i.e. if the left snarl has already extended it's end boundary, we need that information.
                        left_snarl_end_node_id =  max(self.snarl_to_anchors_dictionary[left_snarl_id][0][0].id, self.snarl_to_anchors_dictionary[left_snarl_id][0][-1].id)
                        if settings.DEBUG:
                            print(f"Snarl on the left is {left_snarl_id}, and it's end node is {left_snarl_end_node_id}", flush=True, file=stderr)
                    
                    left_degree = self.graph.get_degree(current_snarl_start_handle, go_left_bool)
                    if settings.DEBUG:
                        print(f"..left_degree of {current_snarl_start_id} is {left_degree}", flush=True, file=stderr)
                    if (
                        left_degree == 2
                        and left_snarl_end_node_id != -1
                        and left_snarl_end_node_id == current_snarl_start_id
                        and self.snarl_to_anchors_dictionary[left_snarl_id][0].bp_occupied_end_node + self.snarl_to_anchors_dictionary[current_snarl_id][0].bp_occupied_start_node == self.graph.get_length(current_snarl_start_handle)
                    ):
                        # try merging to one direction
                        if settings.DEBUG:
                            print(f"    ..Trying _extending_anchors_by_merging in left", flush=True, file=stderr)
                        current_snarl_anchors, snarl_ids_list_idx = self._extending_anchors_by_merging(snarl_ids_sorted, snarl_ids_list_idx, current_snarl_id, left_snarl_id, current_snarl_anchors, extend_left=extend_left, anchors_to_discard=anchors_to_remove, snarl_orientation=snarl_orientation, merging_round=merging_round)
                        if settings.DEBUG:
                            print(f"    ..#anchors returned after merging snarls {current_snarl_id} and {left_snarl_id}: ", len(current_snarl_anchors), flush=True, file=stderr)
                            print(f"    ..new snarl id after merging is: {current_snarl_anchors[0].snarl_id}", flush=True, file=stderr)
                        if current_snarl_anchors[0].snarl_id != current_snarl_id:
                            # TODO: Try to truncate valid_anchors[anchor_idx][1] to read[:4]
                            valid_anchors.extend([[anchor_i, anchor_i.bp_matched_reads] for anchor_i in current_snarl_anchors])
                    min_anchor_length_in_snarl = min([anchor.basepairlength for anchor in current_snarl_anchors])
                    if min_anchor_length_in_snarl >= settings.MIN_ANCHOR_LENGTH:
                        break

                    extend_left = not extend_left
                    current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]
                    current_snarl_end_id = max(current_snarl_anchors[0][0].id, current_snarl_anchors[0][-1].id)
                    if settings.DEBUG:
                        print(f"Current snarl end node is {current_snarl_end_id}", flush=True, file=stderr)
                    current_snarl_end_handle = self.graph.get_handle(current_snarl_end_id)
                    
                    right_snarl_id = snarl_ids_sorted[snarl_ids_sorted.index(current_snarl_id) + 1] if (snarl_ids_sorted.index(current_snarl_id) + 1 < len(snarl_ids_sorted)) else -1 
                    right_snarl_start_node_id = 100000000000000
                    if (
                        self.snarl_to_anchors_dictionary.get(right_snarl_id) != None
                        and len(self.snarl_to_anchors_dictionary[right_snarl_id]) > 0
                    ):
                        right_snarl_start_node_id = min(self.snarl_to_anchors_dictionary[right_snarl_id][0][0].id, self.snarl_to_anchors_dictionary[right_snarl_id][0][-1].id) 
                        if settings.DEBUG:
                            print(f"Snarl on the right is {right_snarl_id}, and it's start node is {right_snarl_start_node_id}", flush=True, file=stderr)

                    right_degree = self.graph.get_degree(current_snarl_end_handle, not go_left_bool)
                    if settings.DEBUG:
                        print(f"..right_degree of {current_snarl_end_id} is {right_degree}", flush=True, file=stderr)
                    if (
                        right_degree == 2
                        and right_snarl_start_node_id != 100000000000000
                        and right_snarl_start_node_id == current_snarl_end_id
                        and self.snarl_to_anchors_dictionary[right_snarl_id][0].bp_occupied_start_node + self.snarl_to_anchors_dictionary[current_snarl_id][0].bp_occupied_end_node == self.graph.get_length(current_snarl_end_handle)
                    ):
                        # current_snarl_id is fetched from list again, as it might have been updated in left-extension
                        current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]
                        if settings.DEBUG:
                            print(f"    ..Trying _extending_anchors_by_merging in right", flush=True, file=stderr)
                        current_snarl_anchors, snarl_ids_list_idx = self._extending_anchors_by_merging(snarl_ids_sorted, snarl_ids_list_idx, current_snarl_id, right_snarl_id, current_snarl_anchors, extend_left=extend_left, anchors_to_discard=anchors_to_remove, snarl_orientation=snarl_orientation, merging_round=merging_round)
                        if settings.DEBUG:
                            print(f"    ..#anchors returned after merging snarls {current_snarl_id} and {right_snarl_id}: ", len(current_snarl_anchors), flush=True, file=stderr)
                            print(f"    ..new snarl id after merging is: {current_snarl_anchors[0].snarl_id}", flush=True, file=stderr)

                        if current_snarl_anchors[0].snarl_id != current_snarl_id:
                            # TODO: Try to truncate valid_anchors[anchor_idx][1] to read[:4]
                            valid_anchors.extend([[anchor_i, anchor_i.bp_matched_reads] for anchor_i in current_snarl_anchors])
                        if settings.DEBUG:
                            print(f"finished extending snarl {current_snarl_id}", flush=True, file=stderr)
                    min_anchor_length_in_snarl = min([anchor.basepairlength for anchor in current_snarl_anchors])
                current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]

            snarl_ids_list_idx += 1

        # now loop over valid_anchors dict to drop all anchors in anchors_to_remove
        for anchor, reads in valid_anchors:
            if isinstance(anchor.snarl_id, str) and "-" in anchor.snarl_id:
                if settings.DEBUG:
                    print(f"snarl {anchor.snarl_id} before adding", flush=True, file=stderr)
            if f"{anchor!r}" not in anchors_to_remove:
                if isinstance(anchor.snarl_id, str) and "-" in anchor.snarl_id:
                    if settings.DEBUG:
                        print(f"snarl {anchor.snarl_id} after adding", flush=True, file=stderr)
                read_info_for_anchor_to_shasta = [read[:4] for read in anchor.bp_matched_reads]
                valid_anchor_extended.append([anchor, read_info_for_anchor_to_shasta])

        return valid_anchor_extended


    def next_handle_iteratee(self, next_boundary):
        self.next_handle_expand_boundary = next_boundary
        # returning False as there is just 1 node connected when the degree is 1.
        return False
    
    def _prepare_snarl_id_chunks_for_parallel_processing(self) -> list[list]:
        """
        Divide the snarl IDs list into chunks. Also return a list of snarl to anchors dictionary for each chunk.
        """

        # TODO: 
        # 1. Generate a #threads vs. runtime plot.
        # 2. Based on the above analysis, set a minimum chunk size. 
        chunk_size = (len(self.snarl_ids_sorted) + self.threads - 1) // self.threads
        
        chunk_snarl_ids_list = [self.snarl_ids_sorted[i:i + chunk_size] for i in range(0, len(self.snarl_ids_sorted), chunk_size)]  
        # egs. [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        
        # generate a list of snarl_to_anchors_dictionary for each chunk
        # chunk_snarl_to_anchors_dict = [{snarl_id: self.snarl_to_anchors_dictionary[snarl_id] for snarl_id in chunk_snarl_ids_list[chunk_idx]} for chunk_idx in range(len(chunk_snarl_ids_list))]  
        # egs. [{1: [anchor1, anchor2], 2: [anchor3, anchor4], 3: [anchor5, anchor6]}, {4: [anchor7, anchor8], 5: [anchor9, anchor10], 6: [anchor11, anchor12]}]

        return chunk_snarl_ids_list
    

    def merge_reliability_checking_results(self, results, file_paths):
        """
        Merge the results from a worker process into the main AlignAnchor instance.
        """
        
        for result in results:
            if settings.OUTPUT_LOGGING_FILES:
                self.snarl_variant_type_dict.update(result["snarl_variant_type_dict"])
                self.snarl_coverage_dict.update(result["snarl_coverage_dict"])
                self.snarl_allelic_coverage_dict.update(result["snarl_allelic_coverage_dict"])
                self.snarl_common_reads_dict.update(result["snarl_common_reads_dict"])
                self.linked_snarls_dictionary.update(result["linked_snarls_dictionary"])
                # self.linked_snarls_compatibility_dict.update(result["linked_snarls_compatibility_dict"])  # This is not correct as it will overwrite the existing dictionary.
                for snarl_id, linked_snarls in result["linked_snarls_compatibility_dict"].items():
                    if snarl_id not in self.linked_snarls_compatibility_dict:
                        self.linked_snarls_compatibility_dict[snarl_id] = {}
                    self.linked_snarls_compatibility_dict[snarl_id].update(linked_snarls)
                self.snarl_read_partitions_dict.update(result["snarl_read_partitions_dict"])
                self.outputs_for_file.extend(result["outputs_for_file"])

            self.reliable_snarls.extend(result["reliable_snarls"])
        
        
        if settings.OUTPUT_LOGGING_FILES:
            # Update the files
            with open(file_paths[0], "w") as f:
                print("snarl_id\tzygosity\tis_reliable\tlinked_snarls", file=f)
                for output in self.outputs_for_file:
                    print(output, file=f)
            
            # Dump the dictionaries
            with open(file_paths[1], "w") as f:
                json.dump(self.snarl_variant_type_dict, f, indent=4)

            with open(file_paths[2], "w") as f:
                json.dump(self.linked_snarls_compatibility_dict, f, indent=4)
            
            with open(file_paths[3], "w") as f:
                json.dump(self.snarl_coverage_dict, f, indent=4)
            
            with open(file_paths[4], "w") as f:
                json.dump(self.snarl_allelic_coverage_dict, f, indent=4)

            with open(file_paths[5], "w") as f:
                json.dump(self.snarl_common_reads_dict, f, indent=4)

            with open(file_paths[6], "w") as f:
                json.dump(self.snarl_read_partitions_dict, f, indent=4)

        # Valid anchors is of format [[anchor, anchor.bp_matched_reads[:4]], [anchor, anchor.bp_matched_reads[:4]], ...]

    
    def dump_valid_anchors(self, extended_out_file_path, anchor_read_tracking_file_path=None, 
                           independent_anchor_read_tracking_file_path=None, reliable_snarls_out_file_path=None, 
                           snarl_variant_type_out_file_path=None, snarl_compatibility_out_file_path=None, snarl_common_reads_out_file_path=None, 
                           snarl_read_partitions_out_file_path=None, snarl_coverage_out_file_path=None, snarl_allelic_coverage_out_file_path=None,
                           snarl_coverage_extended_out_file_path=None, snarl_allelic_coverage_extended_out_file_path=None) -> list:
        
        """
        It iterates over the anchor dictionary. If it finds an anchor with > READS_DEPTH sequences that align to it,
        it adds the list of reads information to the list of anchors to provide as output in json format.
        If sentinel A1 has 2 sequences (S1 and S2) that align to its anchor A1 
        and two sequences (S3 and S4) that align to its anchor A2,
        the valid_anchors list structure is:
                    [ [ [S1A1], [S2A1] ], [ [S3A2], [S4A2] ], [ [S2B1], [S3B1] ] ]
        anchors:       -------A1-------    -------A2-------    -------B1-------


        Parameters
        ----------
        out_file_path : string
            the file of the path were to dump the json file

        Returns
        -------
        bool
        True if iteration has to continue else False
        """

        valid_anchors = []
        valid_anchors_to_extend = []

        for sentinel in self.anchor_reads_dict:
            for id, reads in enumerate(self.anchor_reads_dict[sentinel]):   # A sentinel could have multiple anchors. Those are interated over by the "id"
                if len(reads) > settings.MIN_ANCHOR_READS:
                    anchor = self.sentinel_to_anchor[sentinel][id]
                    snarl_id = anchor.snarl_id
                    self.snarl_to_anchors_dictionary[snarl_id].append(anchor)    # stores snarl to anchors mapping for anchor extension
                    for read in reads:
                        read_id = read[0]
                        if read_id not in self.read_to_snarl_dictionary:
                            self.read_to_snarl_dictionary[read_id] = []
                        self.read_to_snarl_dictionary[read_id].append(anchor.snarl_id)  # stores read IDs and the snarls it passes through (read journey)

        # ## Sort the snarl IDs based on the anchor precedence
        # def anchor_custom_comparator_wrapper(snarl_id1, snarl_id2):
        #     anchor1 = self.snarl_to_anchors_dictionary[snarl_id1][0]
        #     anchor2 = self.snarl_to_anchors_dictionary[snarl_id2][0]
        #     return anchor1.is_preceding_anchor(anchor2)

        # self.snarl_ids_sorted = sorted(list(self.snarl_to_anchors_dictionary.keys()), key=cmp_to_key(anchor_custom_comparator_wrapper))
        # self.snarl_ids_sorted = sorted(list(self.snarl_to_anchors_dictionary.keys()))

        # NOTE: no need to sort the snarls here. Will do sorting on the reliable snarl list later.
        self.snarl_ids_sorted = list(self.snarl_to_anchors_dictionary.keys())
        
        ########### PARALLELIZED: FINDING RELIABLE SNARLS ###########
        t_0 = time.time()
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(f"Processing snarl IDs in parallel with {self.threads} threads...", flush=True, file=stderr)

        # CRITICAL FIX: Remove the C++ PackedGraph object before forking
        # Problem: C++ objects (PackedGraph from bdsg) don't support copy-on-write
        #          When forking, they're copied immediately (~1.6GB per worker)
        # Solution: Remove graph before fork since snarl processing doesn't need it
        # Impact: Reduces memory by ~18 MB per worker
        graph_backup = self.graph
        self.graph = None

        # Set global variable before forking to leverage copy-on-write (avoids pickling)
        global shared_align_anchor
        shared_align_anchor = self

        # Divide the snarl IDs list into chunks
        list_of_chunked_snarl_ids = self._prepare_snarl_id_chunks_for_parallel_processing()
        with multiprocessing.Pool(processes=self.threads, initializer=init_worker_snarl) as pool:
            results = pool.map(process_each_snarl_chunk_in_worker, list_of_chunked_snarl_ids)
        
        # Restore the graph after multiprocessing completes
        self.graph = graph_backup
        
        if settings.DEBUG:
            print("Merging results from worker processes...", flush=True, file=stderr)
        
        if settings.OUTPUT_LOGGING_FILES:
            file_paths = [
                reliable_snarls_out_file_path, 
                snarl_variant_type_out_file_path, 
                snarl_compatibility_out_file_path, 
                snarl_coverage_out_file_path, 
                snarl_allelic_coverage_out_file_path, 
                snarl_common_reads_out_file_path, 
                snarl_read_partitions_out_file_path
            ]
        else:
            file_paths = []
        self.merge_reliability_checking_results(results, file_paths)

        # NOTE:
        # Changelog: Earlier, valid_anchors_from_reliable_snarls was being returned from the merge_reliability_checking_results(...).
        # But our goal is for anchors in valid_anchors_from_reliable_snarls and self.snarl_to_anchors_dictionary[snarl_id] to point to the same underlying anchor objects in memory.
        # So, we recreate the valid_anchors_from_reliable_snarls afresh, using the snarl_ids in self.reliable_snarls, and the self.snarl_to_anchors_dictionary[snarl_id]
        valid_anchors_from_reliable_snarls = []
        self.snarl_ids_sorted = sorted(self.reliable_snarls)

        for snarl_id in self.snarl_ids_sorted:
            for anchor in self.snarl_to_anchors_dictionary[snarl_id]:
                valid_anchors_from_reliable_snarls.append([anchor, [read[:4] for read in anchor.bp_matched_reads]])

        self.runtime_logs.update({"threads": self.threads, "time_for_reliable_snarls_finding": time.time() - t_0})
        
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(f".. Found reliable snarls in {time.time() - t_0}s", flush=True, file=stderr)

        ########### NOT PARALLELIZED ###########
        if settings.DEBUG:
            print(f"######### EXTENDING AND MERGING SNARLS #########", flush=True, file=stderr)
        t_0 = time.time()
        self.valid_anchors_extended = self.extend_and_merge_snarls(valid_anchors=valid_anchors_from_reliable_snarls)       # make sure that it returns serialized anchor object
        
        if settings.DEBUG:
            print(f"######### DUMPING OUTPUTS #########", flush=True, file=stderr)
        
        # Always filter out anchors with less than MIN_ANCHOR_LENGTH and then dump to jsonl
        # FIXME: Make more efficient.
        dump_to_jsonl([[f"{anchor!r}", reads] for anchor, reads in self.valid_anchors_extended if anchor.basepairlength >= settings.MIN_ANCHOR_LENGTH], extended_out_file_path)   # also dumping valid_anchors_extended
        
        
        if settings.OUTPUT_LOGGING_FILES:
            # Populate the snarl_coverage_dict and snarl_allelic_coverage_dict for the extended snarls
            for anchor, reads in self.valid_anchors_extended:
                self.extended_snarl_coverage_dict[anchor.snarl_id] = len(anchor.bp_matched_reads)
                self.extended_snarl_allelic_coverage_dict[anchor.snarl_id] = {
                    idx: len(anchor.bp_matched_reads)
                    for idx, anchor in enumerate(self.snarl_to_anchors_dictionary[anchor.snarl_id])
                }
            
            with open(snarl_coverage_extended_out_file_path, "w") as f:
                json.dump(self.extended_snarl_coverage_dict, f, indent=4)
            
            with open(snarl_allelic_coverage_extended_out_file_path, "w") as f:
                json.dump(self.extended_snarl_allelic_coverage_dict, f, indent=4)
        
            dump_to_jsonl(self.anchor_read_tracking_dict, anchor_read_tracking_file_path)                                 # currently, read drop during snarl merging is not being tracked
            dump_to_jsonl(self.independent_anchor_extension_tracking_dict, independent_anchor_read_tracking_file_path)    # dumping independent anchor extension tracking
        
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(f"Extending and merging snarls took {time.time() - t_0} seconds", flush=True, file=stderr)

        return


    def _find_linked_snarls_for_current_snarl(self, current_snarl_id: str, snarl_list: list, local_snarl_coverage_dict: dict=None, local_snarl_allelic_coverage_dict: dict=None) -> dict:
        """
        Find snarls linked to the current snarl and count the common reads. 
        For S an informative linked snarl T is one such that there exist at least k shared reads
        and in each of S and T the shared reads are partitioned into at least two alleles/anchors.

        Optimization: For checking snarl linkage and compatibility, only check snarls linked by a read 
        """
        linked_snarl_counts = {}

        # Precompute read sets for each anchor in the current snarl
        current_snarl_anchor_sets = [
            {read[0] for read in anchor.bp_matched_reads}
            for anchor in self.snarl_to_anchors_dictionary[current_snarl_id]
        ]
        all_current_reads = set().union(*current_snarl_anchor_sets)

        # Populate the snarl_coverage dictionary
        if local_snarl_coverage_dict is not None:
            local_snarl_coverage_dict[current_snarl_id] = len(all_current_reads)

        # Populate the snarl_allelic_coverage dictionary
        if local_snarl_allelic_coverage_dict is not None:
            local_snarl_allelic_coverage_dict[current_snarl_id] = {
            idx: len(anchor.bp_matched_reads)
            for idx, anchor in enumerate(self.snarl_to_anchors_dictionary[current_snarl_id])
        }

        ## Find snarls linked by current snarl's reads
        potentially_linked_snarls = set(
            snarl_id
            for read_id in all_current_reads
            for snarl_id in self.read_to_snarl_dictionary[read_id]
        )
        
        if settings.DEBUG:
            print(f"For reliability, checking linkage of {current_snarl_id} with {len(potentially_linked_snarls)} snarls", flush=True, file=stderr)

        for other_snarl_id in potentially_linked_snarls:
            if other_snarl_id == current_snarl_id:
                continue  # skip self-comparison

            other_snarl_anchors = self.snarl_to_anchors_dictionary[other_snarl_id]
            # Precompute read sets for other snarl anchors
            other_snarl_anchor_sets = [
                {read[0] for read in anchor.bp_matched_reads}
                for anchor in other_snarl_anchors
            ]
            all_other_reads = set().union(*other_snarl_anchor_sets)

            shared_reads = all_current_reads & all_other_reads
            total_common_reads = len(shared_reads)

            if total_common_reads < settings.MIN_SNARL_LINKAGE_THRESHOLD:
                continue
            
            # Check that shared reads are partitioned into at least two alleles in the current snarl
            current_snarl_partitions = sum(1 for anchor_read_set in current_snarl_anchor_sets if anchor_read_set & shared_reads)
            if current_snarl_partitions < 2:
                continue

            # Check that shared reads are partitioned into at least two alleles in the other snarl
            other_snarl_partitions = sum(1 for anchor_read_set in other_snarl_anchor_sets if anchor_read_set & shared_reads)
            if other_snarl_partitions < 2:
                continue

            # print(f"Found {total_common_reads} common reads between {current_snarl_id} and {other_snarl_id}")
            linked_snarl_counts[other_snarl_id] = total_common_reads

        if settings.DEBUG:
            print(f"Found {len(linked_snarl_counts)} linked snarls for {current_snarl_id}", flush=True, file=stderr)
        
        return linked_snarl_counts


    def _is_other_superset_of_primary(self, primary_sets: list, other_sets: list) -> bool:
        """
        Check if the primary sets are a superset of the other sets.
        """
        for primary_set in primary_sets:
            found = False
            for other_set in other_sets:
                if primary_set.issubset(other_set):
                    found = True
                    break
            if not found:
                return False
        return True


    def _are_unequal_number_of_sets_compatible(self, primary_sets: list, other_sets: list) -> bool:
        """
        Check if two sets are compatible when the cardinality of sets is different.
        """
        # Check if the number of sets is different by more than 1
        if self._is_other_superset_of_primary(primary_sets, other_sets) or self._is_other_superset_of_primary(other_sets, primary_sets):
            return True
        return False
       

    def _are_snarls_compatible(self, primary_snarl: str, other_snarl: str, snarl_read_partitions_dict: dict=None) -> bool:
        """
        Check if two snarls are compatible: 
        S and T linked snarls are consistent if the partition of the shared reads is the "same" in both.
        This means the shared reads should be partitioned identically across the anchors in both snarls.
        For example:
            * Primary snarl has read partitions {1,2,3} and {4,5,6} and other snarl has read partitions {4,5,6} and {1,2,3}, then they are compatible.
            * Primary snarl has read partitions {1,2,3} and {4,5,6} and other snarl has read partitions {4,5,6} and {1,2,3,7}, then they are not compatible.
            * Primary snarl has read partitions {1,2,3}, {4,5,6}, {7,8,9} and other snarl has read partitions {1,2,3,4,5,6}, {7,8,9} then they are not compatible.
        """

        # Collect all read IDs in other_snarl
        other_snarl_reads = {
            read[settings.READ_ID]
            for anchor in self.snarl_to_anchors_dictionary[other_snarl]
            for read in anchor.bp_matched_reads
        }

        # Find common reads between primary and other snarls
        common_reads = {
            read[settings.READ_ID]
            for anchor in self.snarl_to_anchors_dictionary[primary_snarl]
            for read in anchor.bp_matched_reads
            if read[settings.READ_ID] in other_snarl_reads
        }
        # print(f".. {len(common_reads)} Common reads: {common_reads}")

        # Filter both snarls' anchors to include only common reads
        primary_sets = [
            {read[settings.READ_ID] for read in anchor.bp_matched_reads if read[settings.READ_ID] in common_reads}
            for anchor in self.snarl_to_anchors_dictionary[primary_snarl]
        ]
        other_sets = [
            {read[settings.READ_ID] for read in anchor.bp_matched_reads if read[settings.READ_ID] in common_reads}
            for anchor in self.snarl_to_anchors_dictionary[other_snarl]
        ]

        # Remove empty sets (anchors with no common reads)
        primary_sets = [s for s in primary_sets if s]
        other_sets = [s for s in other_sets if s]

        # print(f"Current primary snarl: {primary_snarl}, other snarl: {other_snarl}")
        # print(f"..Primary sets: {primary_sets}")
        # print(f"..Other sets: {other_sets}")

        # Store the partitions for debugging and analysis
        if snarl_read_partitions_dict is not None:
            if (int(primary_snarl.split("-")[0]) if isinstance(primary_snarl, str) else primary_snarl) < (int(other_snarl.split("-")[0]) if isinstance(other_snarl, str) else other_snarl):
                if primary_snarl not in snarl_read_partitions_dict:
                    snarl_read_partitions_dict[primary_snarl] = {}
                if other_snarl not in snarl_read_partitions_dict[primary_snarl]:
                    snarl_read_partitions_dict[primary_snarl][other_snarl] = {
                        "primary": [list(s) for s in primary_sets],
                        "other": [list(s) for s in other_sets]
                    }

        def _are_sets_equal_gtest(primary_sets, other_sets):
            """
            Check if the two sets are permutation equivalent using the G-test.
            """
            if len(primary_sets) != len(other_sets):
                return (False, "False_setsUnequal")

            tangle_matrix = [[len(primary_set & other_set) for other_set in other_sets] for primary_set in primary_sets]
            gtest = GTest(tangle_matrix, settings.DETANGLE_GTEST_EPSILON)
            if not gtest.success or len(gtest.hypotheses) == 0:   # will happen if the tangle matrix is too large (more than 16 entries) or if there are no hypotheses (can only happen when tangle matrix has 0 entries)
                return (False, "False_gtestFailed")
            if not (gtest.hypotheses[0].isForwardInjective() and gtest.hypotheses[0].isBackwardInjective()):    # means that the best hypothes is bijective (both injective and surjective)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                return (False, "False_bestHypothesisNotABijective")
            if gtest.hypotheses[0].G > settings.DETANGLE_MAX_LOG_P:
                return (False, f"False_bestHypothesisGTooHigh {round(gtest.hypotheses[0].G, 2)} > {settings.DETANGLE_MAX_LOG_P}")
            if (len(gtest.hypotheses) > 1) and (gtest.hypotheses[1].G - gtest.hypotheses[0].G < settings.DETANGLE_MIN_LOG_P_DELTA):
                return (False, f"False_hypothesesNotWellSeparated {round(gtest.hypotheses[1].G - gtest.hypotheses[0].G, 2)} < {settings.DETANGLE_MIN_LOG_P_DELTA}")
            
            # # We need to understand why the snarls were compatible. Whether it was exactly [[0,16],[18,0]], i.e. tangle matrix with 0s, or it had some errors, e.g. [[16,2],[0,16]].
            # best_hypothesis = gtest.hypotheses[0]
            # # get indices of False in the best_hypothesis.connectivityMatrix
            # false_indices = [(i, j) for i, row in enumerate(best_hypothesis.connectivityMatrix) for j, value in enumerate(row) if not value]
            # if len(false_indices) == 0:
            #     return (True, "True_tangleMatrixWith0s")
            # else:
            #     return (True, "True_tangleMatrixWithErrors")

            return (True, "True")
        
        def _are_sets_equal_with_error_tolerance(primary_sets, other_sets, error_tolerance=0.1):
            """
            Check if two sets are equal with error tolerance.
            """
            if len(primary_sets) != len(other_sets):
                return (False, "False_setsUnequal")
                # return (self._are_unequal_number_of_sets_compatible(primary_sets, other_sets) if settings.ENABLE_UNEQUAL_SET_COMPATIBILITY else False)
            other_sets_copy = copy.deepcopy(other_sets)
            for primary_set in primary_sets:
                best_matched_intersection_set_size = 0
                best_matched_other_set = None
                for other_set in other_sets_copy:
                    intersection_set = primary_set & other_set 
                    tolerated_error_read_count_primary = len(primary_set) * error_tolerance
                    tolerated_error_read_count_other = len(other_set) * error_tolerance
                    if (len(primary_set) - len(intersection_set) <= tolerated_error_read_count_primary) and (len(other_set) - len(intersection_set) <= tolerated_error_read_count_other):
                        if len(intersection_set) > best_matched_intersection_set_size:
                            best_matched_intersection_set_size = len(intersection_set)
                            best_matched_other_set = other_set
                if best_matched_intersection_set_size == 0:
                    return (False, "False")
                other_sets_copy.remove(best_matched_other_set)
            # Check if all read sets are above a coverage threshold
            for primary_set in primary_sets:
                if len(primary_set) < settings.MIN_READS_FOR_PARTITION_COMPATIBILITY:
                    return (False, "False_lowCov") # if the primary set has less than MIN_READS_FOR_PARTITION_COMPATIBILITY, then the partitions are not compatible
            return (True, "True")

        if settings.USE_GTEST_FOR_PARTITION_COMPATIBILITY:
            is_compatible, desc = _are_sets_equal_gtest(primary_sets, other_sets)
            if is_compatible:
                return (True, "True")
            else:
                return (False, desc)
        
        is_compatible, desc = _are_sets_equal_with_error_tolerance(primary_sets, other_sets, error_tolerance=settings.ERROR_TOLERANCE_IN_COMPATIBILITY_CHECK)
        if is_compatible:
            return (True, "True")
        else:
            return (False, desc)


    def find_reliable_snarls(self, valid_anchors: list, snarl_list: list) -> dict:
        """
        Finds reliable snarls by checking if the current snarl is compatible with 
        >= RELIABLE_SNARL_FRACTION_THRESHOLD of its linked snarls.
        """

        if settings.OUTPUT_LOGGING_FILES:
            local_snarl_variant_type_dict = {}
            local_snarl_coverage_dict = {}
            local_snarl_allelic_coverage_dict = {}
            local_linked_snarls_dictionary = {}
            local_snarl_common_reads_dict = {}
            local_snarl_read_partitions_dict = {}
            outputs_for_file = []

        local_linked_snarls_compatibility_dict = {}
        local_reliable_snarls = []
        # local_valid_anchors_from_reliable_snarls = []
        
        for snarl_id in snarl_list:
            if settings.OUTPUT_LOGGING_FILES:
                # Populate the snarl_variant_type dictionary (SNP or INDEL)
                anchor_sentinel_lengths = []
                for anchor in self.snarl_to_anchors_dictionary[snarl_id]:
                    anchor_sentinel_lengths.append(anchor.sentinel_length)
                if len(set(anchor_sentinel_lengths)) == 1:
                    if anchor_sentinel_lengths[0] == 1:
                        local_snarl_variant_type_dict[snarl_id] = "SNP"
                    else:
                        local_snarl_variant_type_dict[snarl_id] = "MNP"
                else:
                    local_snarl_variant_type_dict[snarl_id] = "INDEL"
            
            #### 1. Find linked snarls and their common read counts
            if settings.OUTPUT_LOGGING_FILES:
                kwargs = {
                        "local_snarl_coverage_dict": local_snarl_coverage_dict,
                        "local_snarl_allelic_coverage_dict": local_snarl_allelic_coverage_dict
                    }
            else:
                kwargs = {}
            linked_snarls_with_counts = self._find_linked_snarls_for_current_snarl(snarl_id, snarl_list, **kwargs)
            
            if settings.OUTPUT_LOGGING_FILES:
                local_snarl_common_reads_dict[snarl_id] = linked_snarls_with_counts
            
            # TODO: Apparantly sorting is needed to make sure the reliable plot doesn't mess up. Handle that since sorting takes time
            linked_snarls_for_current_snarl = sorted(list(linked_snarls_with_counts.keys()))

            if settings.OUTPUT_LOGGING_FILES:
                local_linked_snarls_dictionary[snarl_id] = linked_snarls_for_current_snarl


            if snarl_id not in local_linked_snarls_compatibility_dict:
                local_linked_snarls_compatibility_dict[snarl_id] = {}

            #### 2. Find compatible and incompatible linked snarls for the current snarl
            for linked_snarl_id in linked_snarls_for_current_snarl:
                if linked_snarl_id not in local_linked_snarls_compatibility_dict:
                    local_linked_snarls_compatibility_dict[linked_snarl_id] = {}

                # 2.1. Check if the snarls are compatible
                if settings.OUTPUT_LOGGING_FILES:
                    kwargs = {
                        "snarl_read_partitions_dict": local_snarl_read_partitions_dict
                    }
                else:
                    kwargs = {}
                is_compatible, desc = self._are_snarls_compatible(primary_snarl = snarl_id, other_snarl = linked_snarl_id, **kwargs)
                if is_compatible:
                    local_linked_snarls_compatibility_dict[snarl_id][linked_snarl_id] = True
                    local_linked_snarls_compatibility_dict[linked_snarl_id][snarl_id] = True
                else:
                    if desc != "False":
                        local_linked_snarls_compatibility_dict[snarl_id][linked_snarl_id] = desc
                        local_linked_snarls_compatibility_dict[linked_snarl_id][snarl_id] = desc
                    else:
                        local_linked_snarls_compatibility_dict[snarl_id][linked_snarl_id] = False
                        local_linked_snarls_compatibility_dict[linked_snarl_id][snarl_id] = False

        #### 3. Find whether the current snarl is "reliable" using number of compatilible linked snarls        
        for idx in range(len(snarl_list)):
            snarl_id = snarl_list[idx]
            zygosity = len(self.snarl_to_anchors_dictionary[snarl_id])
            num_compatible_linked_snarls = sum([ 1 for i in local_linked_snarls_compatibility_dict[snarl_id].values() if i == True ])    # calculating compatible linked snarls
            num_non_hom_total_linked_snarls = sum([ 1 for i in local_linked_snarls_compatibility_dict[snarl_id].values()])   # calculating total linked snarls
            fraction_compatible_linked_snarls = (num_compatible_linked_snarls / num_non_hom_total_linked_snarls) if num_non_hom_total_linked_snarls > 0 else 0
            is_reliable = fraction_compatible_linked_snarls > settings.RELIABLE_SNARL_FRACTION_THRESHOLD
            if is_reliable or (zygosity == 1 if settings.ADD_BACK_HOMO_SNARLS else False):
                local_reliable_snarls.append(snarl_id)
            
            if settings.OUTPUT_LOGGING_FILES:
                outputs_for_file.append(f"{snarl_id}\t{zygosity}\t{is_reliable}\t{local_linked_snarls_dictionary[snarl_id]}")

        # local_valid_anchors_from_reliable_snarls = [ele for ele in valid_anchors if ele[0].snarl_id in local_reliable_snarls]

        if settings.OUTPUT_LOGGING_FILES:
            return {
                "snarl_variant_type_dict": local_snarl_variant_type_dict,
                "snarl_coverage_dict": local_snarl_coverage_dict,
                "snarl_allelic_coverage_dict": local_snarl_allelic_coverage_dict,
                "snarl_common_reads_dict": local_snarl_common_reads_dict,
                "linked_snarls_dictionary": local_linked_snarls_dictionary,
                "linked_snarls_compatibility_dict": local_linked_snarls_compatibility_dict,
                "snarl_read_partitions_dict": local_snarl_read_partitions_dict,
                "reliable_snarls": local_reliable_snarls,
                # "valid_anchors_from_reliable_snarls": local_valid_anchors_from_reliable_snarls,
                "outputs_for_file": outputs_for_file
            }
        else:
            return {
                "reliable_snarls": local_reliable_snarls
                # "valid_anchors_from_reliable_snarls": local_valid_anchors_from_reliable_snarls,
            }


    def print_extended_anchor_info(self, out_f) -> None:
        with open(out_f, "w") as f:
            print(f"Sentinel_node\tsnarl_id\tAnchor_length\tAnchor_pos_in_ref_path\tAnchor_path\tAnchor_nodes_copypaste_bandage\tPaths_associated_with_anchor\tbp_matched_reads",file=f)
            # for anchor, _ in self.valid_anchors_extended:
            for anchor, _ in self.valid_anchors_extended:
                print(
                    f"{anchor.get_sentinel_id()}\t{anchor.snarl_id}\t{anchor.basepairlength}\t{anchor.genomic_position}\t{anchor!r}\t{anchor.bandage_representation()}\t{anchor.get_reference_paths()}\t{len([x[0] for x in anchor.bp_matched_reads])}",
                    file=f,
                )


    def print_sentinels_for_bandage(self, file) -> None:
        with open(file, "w") as out_f:
            print("Node,color", file=out_f)
            for anchor, _ in self.valid_anchors_extended:
                for node in anchor:
                    print(f"{node.id},#e25759", file=out_f)


    def dump_snarls_and_anchors_in_reads_dict(self, out_file_path: str) -> None:
        """
        Builds a nested dictionary mapping read IDs to their snarls and anchors,
        then dumps it to a JSONL file.
        """
        snarls_anchors_in_reads_dict = {
            read_id: {
                snarl_id: [
                    f"{anchor!r}"
                    for anchors in self.snarl_to_anchors_dictionary[snarl_id]
                    for anchor in anchors
                ]
                for snarl_id in snarls
            }
            for read_id, snarls in self.read_to_snarl_dictionary.items()
        }

        dump_to_jsonl(snarls_anchors_in_reads_dict, out_file_path)

    
    def dump_dictionary_with_reads_counts(self, out_file_path: str) -> None:
        """
        It writes the anchor dictionary with the count of alinged reads for each anchor

        Parameters
        ----------
        out_file_path : string
            The path to the pkl object that will the dictionary.
        """
        with open(out_file_path, "wb") as out_f:
            pickle.dump(self.sentinel_to_anchor, out_f)


    def processGafLine(self, alignment_l: list, debug_file: str = None):
        """
        It processes an alignment list (the result of parsing an alignment line) to find anchors in the read associated with the alignment. It returns the results to be collected by the caller.
        It:
        1 - walks on the nodes of the path where the read aligns
        2 - if the node where it is standing is a sentinel, checks if there is an anchor that matches the path around the node.
        3 - if it finds one, it then verifies, using the cs tag if the sequence of the reads aligns perfectly (total match) with the portion of the path that is the anchor.
        4 - if so, it captures the read information for that anchor and stops looking for other anchors for this read.
        5 - keeps walking until the aligned path ends.
        """
        # This dictionary will store the results for the given alignment.
        results = {
            "bp_matched_reads": {},
            "anchor_reads": {}
        }

        read_id = alignment_l[settings.READ_POSITION]
        
        walked_length = 0
        if settings.DEBUG:
            print(f"Processing read {read_id}.....", flush=True, file=stderr)

        for position, node_id in enumerate(alignment_l[settings.NODE_POSITION]):

            # Verifying that the nodes coming from the alingment are in the graph I am using
            if not self.graph.has_node(node_id):
                if settings.DEBUG:
                    print(f"THE NODE {node_id} PRESENT IN THE ALIGNMENT IS NOT IN THE PACKED GRAPH.")
                exit(1)

            node_handle = self.graph.get_handle(node_id)
            length = self.graph.get_length(node_handle)

            anchors = self.sentinel_to_anchor.get(node_id)
            
            if anchors:
                for index, anchor in enumerate(anchors):
                    # print(f"For read {read_id}, checking path concordance for anchor {anchor!r}")

                    # an anchor is a list tuple of a list of node handles
                    # and a counter set to 0 at the beginning

                    # scan backward to check that the alignment corresponds to the anchor.
                    alignment_matches_anchor, walk_start, walk_end, relative_strand, walk_start_for_cs_matching, walk_end_for_cs_matching = (
                        verify_path_concordance(
                            position,
                            node_id,
                            alignment_l[settings.NODE_POSITION],
                            alignment_l[settings.ORIENTATION_POSITION],
                            anchor,
                            walked_length
                        )
                    )
                    
                    if settings.DEBUG:
                        print(f"DEBUG: alignment_matches_anchor: {alignment_matches_anchor}, walk_start: {walk_start}, walk_end: {walk_end}, relative_strand: {relative_strand}, walk_start_for_cs_matching: {walk_start_for_cs_matching}, walk_end_for_cs_matching: {walk_end_for_cs_matching}", flush=True, file=stderr)
                    
                    if alignment_matches_anchor:                        
                        x = (
                            anchor,
                            read_id,
                            walk_start,
                            walk_end,
                            alignment_l[settings.CIGAR_POSITION],
                            alignment_l[settings.START_POSITION],
                            alignment_l[settings.END_POSITION],
                            walk_start_for_cs_matching,
                            walk_end_for_cs_matching,
                            alignment_l[settings.READ_START_POS]
                        )

                        is_aligning, read_start, read_end, match_limit, cs_start_pos, cs_end_pos = (
                            verify_sequence_agreement(*x)
                        )
                        # If paths is correct:
                        # I need to append the read info to the anchor.
                        # I need read start and read end of the anchor and the orientation of the read
                        # if (debug_file):
                            # print(f"{read_id},{repr(anchor)},{alignment_matches_anchor},{is_aligning},{match_limit},{cs_start_pos},{cs_end_pos}", file=debug_file)
                        if is_aligning:
                            # print(f" {anchor!r} bp matched")
                            # self.reads_matching_anchor_sequence += 1                          
                            # if not (alignment_l[STRAND_POSITION]):
                            # TODO: Better relative strand calculation. For first read in anchor, store 0 strand and coordinates. Compute the alignment string.
                            # For next read, if the string is same as previous, then strand = 0, else check if it's reverse complement, then strand = 1. If nothing, then report.                            
                            if not (relative_strand):
                                tmp = read_start
                                read_start = alignment_l[settings.R_LEN_POSITION] - read_end
                                read_end = alignment_l[settings.R_LEN_POSITION] - tmp

                            # strand = 0 if alignment_l[STRAND_POSITION] else 1
                            strand = 0 if relative_strand else 1
                            
                            # Store results keyed by anchor identifier
                            anchor_key = (node_id, index)
                            results["bp_matched_reads"][anchor_key] = [[alignment_l[settings.READ_POSITION], strand, read_start, read_end, match_limit, cs_start_pos, cs_end_pos]]
                            results["anchor_reads"][anchor_key] = [[alignment_l[settings.READ_POSITION], relative_strand, read_start, read_end]]

                            # if node_id in [49638724, 49638725, 49638727] and read_id == "c8cb4810-7d6d-42ea-8680-a0483aaabeb1":
                            #     print(f"DEBUG: anchor {anchor!r}: bp_matched_reads = {results['bp_matched_reads'][anchor_key]}", flush=True, file=stderr)

                            break
            
            # adding to the walked length the one of the node I just passed
            walked_length += length
            
        return results, read_id      # After finding all anchors for a read, return the results


def dump_to_jsonl(object, out_file_path: str):
    """
    It dumps the object to json structure.

    Parameters
    ----------
    valid_anchors : list
        the list containing lists of anchors for each sentinel.
    """
    with open(out_file_path, "w", encoding="utf-8") as f:
        json.dump(object, f, ensure_ascii=False, indent=4)


def verify_path_concordance(
    # self,
    alignment_position: int,
    node_id: tuple,
    alignment_node_id_list: list,
    alignment_orientation_list: list,
    anchor: Anchor,
    walked_length: int
) -> list:
    """
    It verifies that the path around the node where the process_alignment function is standing matches the anchor.
    If so, it returns True and returns how many base pairs before and after the start of the sentinel node the sequence alignment has to be a perfect match to validate the anchor.

    Parameters
    ----------
    alignment_position: int
        The position of the current sentinel node being evaluated by the process_alignment function. Is serves as beginning position to locate the sentinel node in the alignment node list.
    node_id: tuple
        Contains the position in the anchor of the sentinel node and the sentinel node
    concordance_orientation: bool
        True if the orientation of the sentinel node is the same between anchor and aligned path, else False
    alignment_node_id_list: list
        The list of node_id corresponding to the alignment path
    alignment_orientation_list: list
        The list of node orientations correspoding to the nodes in alignment_node_id_list
    anchor: list
        The anchor list
    walked_length: int 
        tot basepairs consumed from the beginning of the alignment

    Returns
    -------
    bool
        True if iteration has to continue else False
    already_walked: int
        The basepairs between the start of the sentinel node and the start of the anchor
    to_walk: int
        The basepairs between the start of the sentinel node and the end of the anchor
    relative_strand: bool
        Orientation of read with respect to the anchor path. True if forward, False if reverse

    """
    # DETERMINING THE POSITION OF THE SENTINEL IN THE ANCHOR PATH
    sentinel_position = next(
                    position
                    for position, node in enumerate(anchor)
                    if node.id == node_id
                )

    # DETERMINING THE ORIENTATION OF THE SENTINEL IN THE ANCHOR PATH
    sentinel_orientation = (
                    True if anchor[sentinel_position].orientation else False
                )

    # DETERMINING IF THE ANCHOR PATH AND THE ALIGNMENT ARE CONCORDANT OR REVERSED
    concordance_orientation = (
                    sentinel_orientation == alignment_orientation_list[alignment_position]
                )
    
    # DETERMINING WHERE IN THE ANCHOR NODES LIST THE SENTINEL IS PLACED. THIS IS USED TO KEEP TRACK OF THE WALKED BASEPARIS
    # The "cut" value tells the function how far from the start of the anchor the sentinel is located.
    sentinel_cut = (
        (len(anchor) - 1 - sentinel_position)
        if not concordance_orientation
        else sentinel_position
    )

    # POSITION OF THE ALIGNMENT AT THE BEGINNING OF THE ANCHOR. IF < 0 OR GREATER THAN ALIGNMENT NODES, EXIT.
    alignment_pos = alignment_position - sentinel_cut
    if alignment_pos < 0 or alignment_pos >= len(alignment_node_id_list):
        return (False, 0, 0, -1, 0, 0)
    
    # INITIALZING A LIST WITH ANCHOR LENGTH TO ZERO. TO KEEP TRACK OF THE BASEPAIRS CONSUMED
    basepairs_consumed_list = [0] * len(anchor)

    # # INITIALIZING A LIST OF SENTINEL NODE LENGTHS TO 0.
    # sentinel_list = anchor.get_sentinels()
    # sentinel_bp_consumed_list = [0] * len(sentinel_list)

    # TO SIMPLIFY OPERATIONS, IF THE ANCHOR IS REVERSED COMPARED TO THE PATH, REVERT THE ANCHOR SO SCANNING IS EASIER
    anchor_concordant = anchor[::-1] if not concordance_orientation else anchor[:]
    
    # POSITION IN SCANNING THE ANCHOR
    anchor_pos = 0
    sentinel_pos = 0

    # BASEPAIR RANGE IN THE ALIGNMENT BETWEEN START AND END OF THE ANCHOR
    # alignment_range = (alignment_pos, alignment_pos + len(anchor))

    # nodes_alignment_range = alignment_node_id_list[alignment_range[0]: alignment_range[1]]
    # orientation_alignment_range = alignment_orientation_list[alignment_range[0]: alignment_range[1]]

    # al_string=""
    # for node,orientation_bool in zip(nodes_alignment_range, orientation_alignment_range):
    #     orientation = ">" if orientation_bool else "<"
    #     al_string += orientation + str(node)

    # SCANNING THE ANCHOR AND ALIGNMENT LIST AT THE SAME TIME. EXIT IF ANY ERROR
    anchor_node_orientations_in_read = []

    while anchor_pos < len(anchor_concordant) and alignment_pos < len(
        alignment_node_id_list
    ):
        # check node_id and concordance is the same
        if (
            alignment_node_id_list[alignment_pos]
            != anchor_concordant[anchor_pos].id
        ) or (
            concordance_orientation
            != (
                alignment_orientation_list[alignment_pos]
                == anchor_concordant[anchor_pos].orientation
            )
        ):
            return (False, 0, 0, -1, 0, 0)

        # Store read orientations w.r.t anchor nodes in path
        anchor_node_orientations_in_read.append(alignment_orientation_list[alignment_pos])
        
        #ADDING THE BASEPAIR LENGTHS
        basepairs_consumed_list[anchor_pos] = anchor_concordant[anchor_pos].length

        # #ADDING SENTINEL NODE LENGTHS if anchor_pos points to one of the sentinel nodes 
        # if anchor_concordant[anchor_pos] in sentinel_list:
        #     sentinel_bp_consumed_list[sentinel_pos] = anchor_concordant[anchor_pos].length
        #     sentinel_pos += 1

        # INCREASING POSITION COUNTER
        anchor_pos += 1
        alignment_pos += 1

    if anchor_pos < len(anchor_concordant):
        # didn't finish walking the entire anchor, probably because of alignment_pos < len(alignment_node_id_list)
        return (
            False,
            0,
            0,
            -1,
            0,
            0
        )
    
    # if read_id in ["d59863b0-5ba6-4c3e-ae32-72413357571e", "6e43d5c4-f768-464d-bb18-4220e9d90f5a"] and node_id in [158329263, 158329269]:
    #     print(f"DEBUG: read_id = {read_id}, node_id = {node_id}, walked_length = {walked_length}", flush=True, file=stderr)

    # if read_id in ["3e438e84-b266-4e8e-8649-2b10a30eac7c"] and node_id in [158329254,158329255,158329257]:
    #     print(f"DEBUG: read_id = {read_id}, node_id = {node_id}, walked_length = {walked_length}", flush=True, file=stderr)

    # if read_id in ["4ba88b40-e6f6-449c-9344-ab3e6caf174d"] and node_id in [158265798,158265800]:
    #     print(f"DEBUG: read_id = {read_id}, node_id = {node_id}, walked_length = {walked_length}", flush=True, file=stderr)

    # COMPUTING START AND END OF WALK FOR BASEPAIR SEQUENCE AGREEMENT
    start_walk = walked_length - sum(basepairs_consumed_list[0:sentinel_cut]) + basepairs_consumed_list[0] - (0 if (basepairs_consumed_list[0] == 1) else 1)
    end_walk = walked_length + sum(basepairs_consumed_list[sentinel_cut:]) - basepairs_consumed_list[-1] + (0 if (basepairs_consumed_list[-1] == 1) else 1)
    start_walk_for_cs_matching = start_walk - 1
    end_walk_for_cs_matching = end_walk + 1

    # if read_id in ["4ba88b40-e6f6-449c-9344-ab3e6caf174d"] and node_id in [158265798,158265800]:
    #     print(f"DEBUG: read_id = {read_id}, node_id = {node_id}, start_walk = {start_walk}, end_walk = {end_walk}", flush=True, file=stderr)

    # COMPUTING READ RELATIVE STRAND
    # Simply use concordance_orientation without caching to avoid modifying shared anchor objects
    # (which would trigger copy-on-write in multiprocessing workers)
    relative_strand = concordance_orientation

    return (True, start_walk, end_walk, relative_strand, start_walk_for_cs_matching, end_walk_for_cs_matching)


def verify_sequence_agreement(
    # self,
    anchor: Anchor,
    read_id: str,
    anchor_bp_start: int,
    anchor_bp_end: int,
    cs_walk: list,
    start_in_path: int,
    end_in_path: int,
    walk_start_for_cs_matching: int,
    walk_end_for_cs_matching: int,
    intialise_walked_in_the_sequence_to: int
):
    """
    It uses the parsed cs tag from the gaf to verify that the anchor and the path match at the sequence level.

    Parameters
    ----------
    anchor_bp_start: int
        Start position of the anchor in the path
    anchor_bp_end: int
        End position of the anchor in the path
    cs_walk: list
        the cs tag operations structured as steps in a list
    start_in_path: int
        The alingment start in the path (from gaf)
    end_in_path: int
        The alingment end in the path (from gaf)
    Returns
    -------
    boool
        True if the path and read sequences matches in the anchor section
    walked_in_the_sequence - diff_start: int
        The start of the anchor in the read / 0 if does not match completely
    walked_in_the_sequence - diff_end: int
        The end of the anchor in the read / 0 if does not match completely
    """
    
    print_to_debug = False

    # If anchor overflows the alingment, it is not valid
    if anchor_bp_end > end_in_path or anchor_bp_start < start_in_path or anchor_bp_end < anchor_bp_start:
        return (False, 0, 0, 0, 0, 0)

    walked_in_the_sequence: int = (
        intialise_walked_in_the_sequence_to  # I need this to keep track of anchor position in the sequence
    )
    walked_in_the_path: int = (
        start_in_path  # I need this to keep track of my walk in the path
    )
    allow_seq_diff: bool = (
        True  # I need this to control no variation between anchor and sequence is present. Starting with True, setting to False when walking on anchor coordinates
    )
    if print_to_debug:
        print(f"DEBUG: Initially, walked_in_the_path = {walked_in_the_path}, walked_in_the_sequence = {walked_in_the_sequence}, allow_seq_diff = {allow_seq_diff}")

    # When walking on alingment. Path length is calculated as 'equal + subst + delition'
    # When walking on alingment. Read length is calculated as 'equal+subst+insertion'
    # For the moment, strand can be assumed as positive
    total_matched_bps = 0
    for step in cs_walk:

        if print_to_debug:
            print(f"DEBUG: walked_in_the_path = {walked_in_the_path}, walked_in_the_sequence = {walked_in_the_sequence}")
            print(f"DEBUG: step = {step}")
        if step[0] == "+":
            walked_in_the_sequence += step[1]
        elif step[0] == ":":
            walked_in_the_sequence += step[1]
            walked_in_the_path += step[1]
        elif step[0] == "-":
            walked_in_the_path += step[1]
        elif step[0] == "*":
            walked_in_the_sequence += step[1]
            walked_in_the_path += step[1]

        if walked_in_the_path > anchor_bp_start and allow_seq_diff:
            # I passed the start of the anchor and I was on a difference step. Anchor not good
            if step[0] != ":":
                return (False, 0, 0, 0, 0, 0)
            if print_to_debug:
                print(f"DEBUG: Just passed start of anchor, walked_in_the_path = {walked_in_the_path}, walked_in_the_sequence = {walked_in_the_sequence}")
            # If I passed on a equal step, it is ok. I set allow_differences to false and go on. But before I check if I have surpassed the end of the anchor. If yes return true.
            total_matched_bps = step[1]
            if print_to_debug:
                print(f"DEBUG: Inside passed start of anchor, total_matched_bps in this cs_step (same as step size) = {total_matched_bps}")
            if walked_in_the_path >= walk_end_for_cs_matching:
                if print_to_debug:
                    print(f"DEBUG: End of anchor in same cs_step as start of anchor, as walked_in_the_path = {walked_in_the_path}, walked_in_the_sequence = {walked_in_the_sequence}, walk_end_for_cs_matching = {walk_end_for_cs_matching}")
                diff_start = walked_in_the_path - anchor_bp_start
                diff_end = walked_in_the_path - anchor_bp_end
                if print_to_debug:
                    print(f"DEBUG: diff_start = {diff_start}, diff_end = {diff_end}")
                if print_to_debug:
                    print(f"DEBUG: returning True, read_start = {walked_in_the_sequence - diff_start}, read_end = {walked_in_the_sequence - diff_end}, match_limit = {total_matched_bps}, cs_left_avail = {total_matched_bps - diff_start}, cs_right_avail = {diff_end}")
                # I add a + 1 in the read_end position because of Shasta requirement that the interval is open at the end. The end id in the sequence is of the first nucleotide after the anchor
                # TODO: Currently end_node_pos causes a gap when anchor end node is even #base-pairs.
                return (
                    True,
                    walked_in_the_sequence - diff_start,
                    walked_in_the_sequence - diff_end,
                    total_matched_bps,
                    total_matched_bps - diff_start,
                    diff_end
                )
            else:
                allow_seq_diff = False  # go to the next step

        # Walking in the anchor section and found a diff
        elif not (allow_seq_diff) and step[0] != ":":
            return (False, 0, 0, 0, 0, 0)

        # Walking inside an anchor (that spans multiple ":" tuples
        elif walked_in_the_path > anchor_bp_start and walked_in_the_path < walk_end_for_cs_matching and step[0] == ":" and (not allow_seq_diff):
            total_matched_bps += step[1]

        # I passed the end of the anchor and there was no difference
        elif walked_in_the_path >= walk_end_for_cs_matching:
            total_matched_bps += step[1]
            diff_start = walked_in_the_path - anchor_bp_start
            diff_end = walked_in_the_path - anchor_bp_end
            # I add a + 1 in the read_end position because of Shasta requirement that the interval is open at the end. The end id in the sequence is of the first nucleotide after the anchor
            return (
                True,
                walked_in_the_sequence - diff_start,
                walked_in_the_sequence - diff_end,
                total_matched_bps,
                total_matched_bps - diff_start,
                diff_end
            )

        elif walked_in_the_path > end_in_path:
            return (False, 0, 0, 0, 0, 0)
    return (False, 0, 0, 0, 0, 0)
