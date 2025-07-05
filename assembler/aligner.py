import json
from sys import stderr
import pickle
import os.path
from sys import exit
from collections import defaultdict
import copy

from bdsg.bdsg import PackedGraph
from assembler.anchor import Anchor
from assembler.node import Node
from assembler.constants import (
    READ_POSITION,
    R_LEN_POSITION,
    STRAND_POSITION,
    START_POSITION,
    END_POSITION,
    NODE_POSITION,
    ORIENTATION_POSITION,
    CIGAR_POSITION,
    MIN_ANCHOR_READS,
    MIN_ANCHOR_LENGTH,
    READ_POSITION,
    HET_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING,
    HOMO_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING,
    MIN_READS_REQUIRED_FOR_MERGING_R0,
    MIN_READS_REQUIRED_FOR_MERGING_R1,
    FRACTION_READS_FOR_SNARL_BOUNDARY_EXTENTION,
    MIN_READS_REQUIRED_FOR_BOUNDARY_EXTENSION,
    READ_ID,
    READ_STRAND,
    ANCHOR_START,
    ANCHOR_END,
    CS_LEFT_AVAIL,
    CS_RIGHT_AVAIL,
    MIN_ANCHOR_READCOV,
    DROP_FRACTION
)
from assembler.anchor_coverage import AnchorCoverage

# read_compare = 'm64015_190920_185703/40699337/ccs' #'m64012_190921_234837/35261139/ccs'#'m64012_190921_234837/96012404/ccs'

class AlignAnchor:

    def __init__(self) -> None:
        # useful initialization objects
        self.graph = PackedGraph()
        self.snarl_to_anchor_reads_dictionary = defaultdict(list)
        self.snarl_to_anchors_dictionary = defaultdict(list)
        self.sentinel_to_anchor: dict = dict()
        self.anchor_reads_dict: dict = dict()
        self.reads_matching_anchor_path: int = 0
        self.reads_matching_anchor_sequence: int = 0
        self.next_handle_expand_boundary = None
        self.node_orientations_in_anchor_for_read1 = []
        self.anchor_coverage = AnchorCoverage()  # Add coverage tracking
        self.anchor_read_tracking_dict = {} # Add coverage tracking for anchors
        

    def build(self, dict_path: str, packed_graph_path: str) -> None:

        # loading dictionary
        with open(dict_path, 'rb') as in_f:
            self.sentinel_to_anchor = pickle.load(in_f)

        #loading packedgraph
        self.graph.deserialize(packed_graph_path)

        # initializing output dictionary
        for sentinel, anchors in self.sentinel_to_anchor.items():
            self.anchor_reads_dict[sentinel] = [[] for _ in range(len(anchors))]

        # for sentinel in self.sentinel_to_anchor:
        #     print(f"S_T_A {sentinel} = {self.sentinel_to_anchor[sentinel]}")
        #     break

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
        anchors_to_discard : list
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

        print("current snarl being extended:- ", current_snarl_id)
        print("surrounding snarls:- ", snarl_ids_sorted_list_up_to_date[snarl_ids_sorted_list_iterator_idx-5:min(len(snarl_ids_sorted_list_up_to_date), snarl_ids_sorted_list_iterator_idx+5)])
        
        # calculating if after merging, enough anchors (>=2) will remain 
        cnt_anchors_with_sufficient_read_overlap = 0
        # other_anchor = self.snarl_to_anchors_dictionary[other_snarl_id][0]

        # When merging in round=1 (i.e. lower MIN_READS_REQUIRED_FOR_MERGING value, hence less confidence), merge only if current anchor is heterozygous. 
        # Else, there is no point in sacrificing read coverage by such an amount, if it doesn't help in phasing. 
        if merging_round == 1 and len(current_snarl_anchors) < 2:
            return (current_snarl_anchors, snarl_ids_sorted_list_iterator_idx)
        MIN_READS_REQUIRED_FOR_MERGING = MIN_READS_REQUIRED_FOR_MERGING_R0 if merging_round == 0 else MIN_READS_REQUIRED_FOR_MERGING_R1
        
        for other_anchor in self.snarl_to_anchors_dictionary[other_snarl_id]:
            for anchor in current_snarl_anchors:
                common_paths = set(anchor.reference_paths_covered).intersection(set(other_anchor.reference_paths_covered))
                read_ids_current_anchor = [read[READ_POSITION] for read in anchor.bp_matched_reads]
                read_ids_other_anchor = [read[READ_POSITION] for read in other_anchor.bp_matched_reads]
                common_reads_ids = set(read_ids_current_anchor).intersection(set(read_ids_other_anchor))
                for read in anchor.bp_matched_reads:
                    if read[READ_POSITION] in common_reads_ids:
                        extra_bps = read[CS_LEFT_AVAIL] if extend_left else read[CS_RIGHT_AVAIL]
                        if extra_bps < 1:
                            common_reads_ids.remove(read[READ_POSITION])
                for read in other_anchor.bp_matched_reads:
                    if read[READ_POSITION] in common_reads_ids:
                        extra_bps = read[CS_RIGHT_AVAIL] if extend_left else read[CS_LEFT_AVAIL]
                        if extra_bps < 1:
                            common_reads_ids.remove(read[READ_POSITION])

                print(f"for anchor {anchor!r} and other anchor {other_anchor!r}, common reads between them are {len(common_reads_ids)}")                 
                if len(common_paths) > 0 and (len(common_reads_ids) > MIN_READS_REQUIRED_FOR_MERGING):    # meaning we can create an anchor with this combination
                    cnt_anchors_with_sufficient_read_overlap += 1
                    if cnt_anchors_with_sufficient_read_overlap > 1:
                        break
                    # extendable_anchors_in_current_snarl.append(anchor)
            if cnt_anchors_with_sufficient_read_overlap > 1:
                break
        
        # [55, 55, 15]
        # l_snarl (hom) = 110
        # r_snarl (hom) = 75
        # l --> [50, 40, 10]
        # r --> [34, 33]

        if cnt_anchors_with_sufficient_read_overlap < 2:
            print(f"Failed to merge snarls {current_snarl_id} and {other_snarl_id} because of insufficient anchors after merging")
            return (current_snarl_anchors, snarl_ids_sorted_list_iterator_idx)
        
        # meaning we found that after extension, heterozygosity of current snarl will be maintained.
        # extend this side
        new_anchors_after_merging = []
        # loop over all possible combinations of anchor_in_current_snarl x anchor_in_other_snarl
        for anchor in current_snarl_anchors:
            for other_anchor in self.snarl_to_anchors_dictionary[other_snarl_id]:
                common_paths = set(anchor.reference_paths_covered).intersection(set(other_anchor.reference_paths_covered))
                read_ids_current_anchor = [read[READ_POSITION] for read in anchor.bp_matched_reads]
                read_ids_other_anchor = [read[READ_POSITION] for read in other_anchor.bp_matched_reads]
                common_reads_ids = set(read_ids_current_anchor).intersection(set(read_ids_other_anchor))
                for read in anchor.bp_matched_reads:
                    if read[READ_POSITION] in common_reads_ids:
                        extra_bps = read[CS_LEFT_AVAIL] if extend_left else read[CS_RIGHT_AVAIL]
                        if extra_bps < 1:
                            print(f"read {read[READ_POSITION]} rejected because of possibility of insertion/deletion between anchors {anchor!r} and {other_anchor!r}.")
                            common_reads_ids.remove(read[READ_POSITION])
                for read in other_anchor.bp_matched_reads:
                    if read[READ_POSITION] in common_reads_ids:
                        extra_bps = read[CS_RIGHT_AVAIL] if extend_left else read[CS_LEFT_AVAIL]    # THIS SHOULD BE CS_RIGHT_AVAIL IF EXTEND_LEFT IS FALSE
                        if extra_bps < 1:
                            print(f"read {read[READ_POSITION]} rejected because of possibility of insertion/deletion between anchors {anchor!r} and {other_anchor!r}.")
                            common_reads_ids.remove(read[READ_POSITION])

                if len(common_paths) > 0 and (len(common_reads_ids) > MIN_READS_REQUIRED_FOR_MERGING):    # meaning we can create an anchor with this combination
                    print(f"merging anchors {anchor!r} (current_anchor) and {other_anchor!r} (other_anchor) in", "left extension." if extend_left==True else "right extension.", end=" ")
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
                        if read[READ_ID] in common_reads_ids:
                            common_bp_matched_reads[read[READ_ID]] = read
                    for read in other_anchor.bp_matched_reads:
                        if read[READ_ID] in common_reads_ids:
                            print(f"processing read {read[READ_ID]}")
                            unpacked_read_id, unpacked_strand, unpacked_start, unpacked_end, unpacked_match_limit, unpacked_cs_left, unpacked_cs_right = common_bp_matched_reads[read[READ_ID]]
                            print(f"current anchor boundary before merge: {anchor!r} : {unpacked_start} - {unpacked_end}")
                            print(f"other anchor boundary before merge: {other_anchor!r} : {read[ANCHOR_START]} - {read[ANCHOR_END]}")
                            unpacked_start = min(read[ANCHOR_START], unpacked_start)
                            unpacked_end = max(read[ANCHOR_END], unpacked_end)
                            unpacked_cs_left = min(read[CS_LEFT_AVAIL], unpacked_cs_left)
                            unpacked_cs_right = min(read[CS_RIGHT_AVAIL], unpacked_cs_right)
                            common_bp_matched_reads[read[READ_ID]] = [unpacked_read_id, unpacked_strand, unpacked_start, unpacked_end, unpacked_match_limit, unpacked_cs_left, unpacked_cs_right]
                            print(f"anchor boundary AFTER merge: {new_anchor!r} : {unpacked_start} - {unpacked_end}")
                    
                    common_bp_matched_reads_list = list(common_bp_matched_reads.values())
                    new_anchor.bp_matched_reads = sorted(common_bp_matched_reads_list, key=lambda read: read[READ_ID])    # set(anchor.bp_matched_reads).intersection(set(other_anchor.bp_matched_reads))
                    new_anchors_after_merging.append(new_anchor)

        # remove all anchors from both current and other snarls, as they are now replaced by new anchors having new snarl name
        anchors_to_discard.extend(current_snarl_anchors)
        anchors_to_discard.extend(self.snarl_to_anchors_dictionary[other_snarl_id])
        # add the new snarl and its anchors in the snarl_to_anchors_dictionary
        new_snarl_id_after_merge = new_anchors_after_merging[0].snarl_id
        print(f"Inside snarl merging, have >=2 new anchors after merging {current_snarl_id} and {other_snarl_id}. New snarl id should be: {new_snarl_id_after_merge}")
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

        print(f"    ...currently checking extension in node ID {self.graph.get_id(node_handle)} which was extend_left = {extend_left} direction.")

        num_anchors_remaining_after_extension = 0
        anchors_to_extend = []
        for anchor in current_snarl_anchors:
            total_bp_matched_reads = len(anchor.bp_matched_reads)
            print(f"    ...anchor pre-extension is: {anchor!r}, total_bp_matched_reads = {total_bp_matched_reads}")
            common_bp_matched_reads = []
            num_reads_that_can_be_extended = 0
            
            bp_added_upon_extention = (self.graph.get_length(current_snarl_boundary_handle) + 1) // 2 + (self.graph.get_length(node_handle_to_extend_to)) // 2 + 1
            # print(f"  bp_added_upon_extention = {bp_added_upon_extention}")

            for read in anchor.bp_matched_reads:
                read_id, read_strand, anchor_start, anchor_end, match_limit, cs_avail_left, cs_avail_right = read
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
            if (num_reads_that_can_be_extended/total_bp_matched_reads >= FRACTION_READS_FOR_SNARL_BOUNDARY_EXTENTION) and (num_reads_that_can_be_extended >= MIN_READS_REQUIRED_FOR_BOUNDARY_EXTENSION):
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
                print(f"    ....INSIDE: New snarl ID of anchor {anchor!r} after extending is {anchor.snarl_id}")
                anchor.bp_matched_reads = common_bp_matched_reads
            
            extended_anchors = [anchor for anchor, _ in anchors_to_extend]
            print(f"    ....INSIDE: extended anchors are {extended_anchors}")
            return extended_anchors
        
        else:
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
                print(f"bp_occupied_next_node = {bp_occupied_next_node}, bp_occupied_other_snarl_boundary_node = {bp_occupied_other_snarl_boundary_node}")
                bp_available_for_extension = self.graph.get_length(next_node_handle) - bp_occupied_next_node - bp_occupied_other_snarl_boundary_node     # bps available in current new node 
                print(f"bp_available_for_extension = {bp_available_for_extension}")
                # TODO: think about base case, when bp_available_for_extension = 0
                #     set a flag here representing that extension can no longer be done after this node
                cant_extend_more = True
            else:    # means this node isn't a boundary node for another snarl, so at max, it can be consumed completely
                print(f"bp_occupied_next_node = {bp_occupied_next_node}")
                bp_available_for_extension = self.graph.get_length(next_node_handle) - bp_occupied_next_node
                print(f"bp_available_for_extension = {bp_available_for_extension}")

            # find overall bps to extend by doing min of bp_available_for_extension with all bps in per_anchor_max_bps_to_extend
            min_bp_among_cs_lines = min(per_anchor_max_bps_to_extend)
            if min_bp_among_cs_lines <= bp_available_for_extension:    # this means cs_avail will be completely consumed in this node 
                cant_extend_more = True
            final_bp_count_added_in_current_iteration = min(min_bp_among_cs_lines, bp_available_for_extension)

            # DOES a_current_snarl_anchor.basepairlength GET CORRECTLY UPDATED IN THIS WHILE LOOP? (yes)
            min_basepairlength_among_snarl_anchors = min(anchor.basepairlength for anchor in current_snarl_anchors)
            if (final_bp_count_added_in_current_iteration + min_basepairlength_among_snarl_anchors > MIN_ANCHOR_LENGTH):
                final_bp_count_added_in_current_iteration = max(0, MIN_ANCHOR_LENGTH - min_basepairlength_among_snarl_anchors)    # THIS SHOULD BE MIN_ANCHOR_LENGTH NOT MIN_ANCHOR_READS
                print(f"We are in extension iteration 2 and final_bp_count_added_in_current_iteration = {final_bp_count_added_in_current_iteration}")
                cant_extend_more = True

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
                    print(f"bp_occupied_start_node of anchor {anchor!r} is {anchor.bp_occupied_start_node}")
                else:
                    anchor.bp_occupied_end_node += final_bp_count_added_in_current_iteration
                    print(f"bp_occupied_end_node of anchor {anchor!r} is {anchor.bp_occupied_end_node}")

                anchor.compute_bp_length()

                # update bp_matched_reads (also update cs_avail_left/right accordingly)
                new_bp_matched_reads = []
                for read in anchor.bp_matched_reads:
                    if extend_left:
                        new_cs_avail_idx = (CS_LEFT_AVAIL if read[READ_STRAND] == 0 else CS_RIGHT_AVAIL)
                        read[ANCHOR_START] = read[ANCHOR_START] - final_bp_count_added_in_current_iteration
                    else:
                        new_cs_avail_idx = (CS_RIGHT_AVAIL if read[READ_STRAND] == 0 else CS_LEFT_AVAIL)
                        read[ANCHOR_END] = read[ANCHOR_END] + final_bp_count_added_in_current_iteration
                    new_cs_avail = read[new_cs_avail_idx] - final_bp_count_added_in_current_iteration
                    if new_cs_avail >= 0:
                        read[new_cs_avail_idx] = new_cs_avail
                        new_bp_matched_reads.append(read)
                    else:
                        if current_snarl_id not in self.anchor_read_tracking_dict:
                            self.anchor_read_tracking_dict[current_snarl_id] = dict()
                        if f"{anchor!r}" not in self.anchor_read_tracking_dict[current_snarl_id]:
                            self.anchor_read_tracking_dict[current_snarl_id][f"{anchor!r}"] = {}
                        if extension_iteration not in self.anchor_read_tracking_dict[current_snarl_id][f"{anchor!r}"]:
                            self.anchor_read_tracking_dict[current_snarl_id][f"{anchor!r}"][extension_iteration] = []
                        # add this read to anchor_read_tracking_dict
                        self.anchor_read_tracking_dict[current_snarl_id][f"{anchor!r}"][extension_iteration].append(read[READ_ID])
                anchor.bp_matched_reads = new_bp_matched_reads

                # update anchor.compute_bp_length() to have correct calculation for boundary nodes
                # if anchor.basepairlength >= MIN_ANCHOR_LENGTH:
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
            allowed_read_drops = int(DROP_FRACTION * current_anchor_readcov)
            if current_anchor_readcov <= MIN_ANCHOR_READCOV:
                return 0
            elif current_anchor_readcov - allowed_read_drops <= MIN_ANCHOR_READCOV:
                return current_anchor_readcov - MIN_ANCHOR_READCOV
            else:
                return allowed_read_drops

        else:
            ##### ITERATION-2: RELAXED READ-DROPS, BUT PUSHES FOR TOUCHING MIN_ANCHOR_LENGTH
            return max(0, current_anchor_readcov - MIN_ANCHOR_READCOV)

            
            # diff_readcov = current_anchor_readcov - MIN_ANCHOR_READCOV

            # if diff_readcov < 0:
            #     return 0
            # elif diff_readcov < MAX_READ_DROPS_ALLOWED:
            #     return diff_readcov
            # else:
            #     return MAX_READ_DROPS_ALLOWED


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
        for anchor in current_snarl_anchors:
            current_snarl_anchor_readcov.append(len(anchor.bp_matched_reads))
        print(f"...current_snarl_anchor_readcov is {current_snarl_anchor_readcov}")

        allowed_read_drop_counts = [0] * len(current_snarl_anchors)    # stores for each anchor, how many reads can be dropped; 0 in the no_drop extension
        if extension_iteration != 0:
            allowed_read_drop_counts = [self._get_max_read_drop(current_anchor_readcov, extension_iteration) for current_anchor_readcov in current_snarl_anchor_readcov]
        allowed_read_drop_counts_iterator = iter(allowed_read_drop_counts)
        per_anchor_max_bps_to_extend_left = []    # storing for each anchor, max cs_avail for extension in the left_direction 
        for anchor in current_snarl_anchors:    # for left extension
            per_read_cs_avail_list = [(read[CS_LEFT_AVAIL] if (read[READ_STRAND] == 0) else read[CS_RIGHT_AVAIL]) for read in anchor.bp_matched_reads]
            current_allowed_read_drop_counts = next(allowed_read_drop_counts_iterator)
            per_anchor_max_bps_to_extend_left.append(self._get_max_cs_avail_in_anchor(per_read_cs_avail_list, current_allowed_read_drop_counts))
        if snarl_ids_list_idx > 0:
            self._try_extension(current_snarl_anchors, current_snarl_id, snarl_ids_sorted[snarl_ids_list_idx - 1], anchors_to_discard, per_anchor_max_bps_to_extend_left, extend_left=True, extension_iteration = extension_iteration)   # for no_drop left extension
        

        current_snarl_anchor_readcov = []    # storing read coverage of each anchor
        for anchor in current_snarl_anchors:
            current_snarl_anchor_readcov.append(len(anchor.bp_matched_reads))
        print(f"...current_snarl_anchor_readcov is {current_snarl_anchor_readcov}")

        allowed_read_drop_counts = [0] * len(current_snarl_anchors)    # stores for each anchor, how many reads can be dropped; 0 in the no_drop extension
        if extension_iteration != 0:
            allowed_read_drop_counts = [self._get_max_read_drop(current_anchor_readcov, extension_iteration) for current_anchor_readcov in current_snarl_anchor_readcov]
        allowed_read_drop_counts_iterator = iter(allowed_read_drop_counts)
        per_anchor_max_bps_to_extend_right = []    # storing for each anchor, max cs_avail for extension in the right_direction 
        for anchor in current_snarl_anchors:    # for right extension
            per_read_cs_avail_list = [(read[CS_RIGHT_AVAIL] if (read[READ_STRAND] == 0) else read[CS_LEFT_AVAIL]) for read in anchor.bp_matched_reads]
            current_allowed_read_drop_counts = next(allowed_read_drop_counts_iterator)
            per_anchor_max_bps_to_extend_right.append(self._get_max_cs_avail_in_anchor(per_read_cs_avail_list, current_allowed_read_drop_counts))
        
        if snarl_ids_list_idx < len(snarl_ids_sorted) - 1:
            self._try_extension(current_snarl_anchors, current_snarl_id, snarl_ids_sorted[snarl_ids_list_idx + 1], anchors_to_discard, per_anchor_max_bps_to_extend_right, extend_left=False, extension_iteration = extension_iteration)   # for no_drop left extension

    
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
            print(f"Processing snarl ID: {current_snarl_id}")
            current_snarl_anchors = self.snarl_to_anchors_dictionary[current_snarl_id]
            print(f"..Running _extending_snarl_boundaries of snarl {current_snarl_id} with {len(current_snarl_anchors)} anchors")
            min_basepairlength_among_snarl_anchors = min([anchor.basepairlength for anchor in current_snarl_anchors])
            if min_basepairlength_among_snarl_anchors >= MIN_ANCHOR_LENGTH:
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
            snarl_ids_list_idx += 1


    def extend_and_merge_snarls(self, valid_anchors: list) -> list:
        """
        This function performs snarl boundary extension and merging
        
        """
        # valid_anchor_extended = []
        # newsnarls_to_anchors_dictionary = {}
        anchors_to_remove = []   # {(snarl_id, anchor)}
        snarl_ids_sorted = sorted(list(self.snarl_to_anchors_dictionary.keys()))
        
        ### First, performing perfect bp match extension (no read drop allowed) for all snarls
        print(f"#### RUNNING EXTENSION WITH NO DROPS AND ALLOWED DROPS FOR HET ANCHORS ONLY ####")
        self._helper_extension_loop(snarl_ids_sorted, anchors_to_remove, extension_iteration=0, is_het_round=True)
        self._helper_extension_loop(snarl_ids_sorted, anchors_to_remove, extension_iteration=1, is_het_round=True)

        print(f"#### RUNNING EXTENSION WITH NO DROPS AND ALLOWED DROPS FOR HOM ANCHORS ONLY ####")
        self._helper_extension_loop(snarl_ids_sorted, anchors_to_remove, extension_iteration=0, is_het_round=False)
        self._helper_extension_loop(snarl_ids_sorted, anchors_to_remove, extension_iteration=1, is_het_round=False)
        
        print(f"#### RUNNING EXTENSION WITH MORE DROPS ALLOWED, FIRST FOR HET, AND THEN HOM ANCHORS ####")
        self._helper_extension_loop(snarl_ids_sorted, anchors_to_remove, extension_iteration=2, is_het_round=True)
        self._helper_extension_loop(snarl_ids_sorted, anchors_to_remove, extension_iteration=2, is_het_round=False)
            

        print(f"#### TRY TO MERGE SHORTER ANCHORS ######")
        valid_anchors = self.merge_anchors(valid_anchors, anchors_to_remove, snarl_ids_sorted, merging_round=0)
        # valid_anchors = self.merge_anchors(valid_anchors, anchors_to_remove, snarl_ids_sorted, merging_round=1)

        print(f"#### MERGING FINISHED... ######")

        for idx in range(len(valid_anchors)):
            anchor = valid_anchors[idx][0]
            reads = [read[:4] for read in anchor.bp_matched_reads]
            valid_anchors[idx][1] = reads
        
        return valid_anchors    #### change this later to calculate valid_anchors_extended, when we will have anchor drops because of merging


    def merge_anchors(self, valid_anchors: list, anchors_to_remove: list, snarl_ids_sorted: list, merging_round: int) -> list:
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
        anchors_to_remove = []   # {(snarl_id, anchor)}
        snarl_orientation = True

        # snarl_ids_sorted = sorted(list(self.snarl_to_anchors_dictionary.keys()))
        snarl_ids_list_idx = 0    # we need to use this index counter (and can't simply use an iterator) because we will be inserting/deleting the snarl_ids_sorted list on the go
        # Note: Remember to update the snarl_ids_list_idx appropriately when merging snarls
        while snarl_ids_list_idx < len(snarl_ids_sorted):
            current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]
            last_snarl_id = -1
            current_snarl_anchors = self.snarl_to_anchors_dictionary[current_snarl_id]
            min_anchor_length_in_snarl = min([anchor.basepairlength for anchor in current_snarl_anchors])
            # if len(current_snarl_anchors) > 1 and min_anchor_length_in_snarl < MIN_ANCHOR_LENGTH:

            while (
                len(current_snarl_anchors) > 1 
                and (min_anchor_length_in_snarl < MIN_ANCHOR_LENGTH
                and last_snarl_id != current_snarl_id)
            ):
                current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]
                last_snarl_id = current_snarl_id

                # recalculating left and right nodes in graph for extension/merging
                if snarl_orientation:
                    # even if 0-th anchor of snarl is reversed, snarl_start and snarl_end should be in increasing order of node ids (if snarl_orientation = True)
                    print(f"########################", end="\n")
                    print(f"Processing SNARL ID: {current_snarl_id}")
                    print(f"####### Trying left 1-degree node extension")
                    current_snarl_start_id = min(current_snarl_anchors[0][0].id, current_snarl_anchors[0][-1].id)   # change this variable to current_snarl_start_node_id
                    print(f"Current snarl start node is {current_snarl_start_id}")
                    go_left_bool = True
                    current_snarl_start_handle = self.graph.get_handle(current_snarl_start_id)
                    extend_left = True

                    # Get left snarl's ID
                    left_snarl_id = snarl_ids_sorted[snarl_ids_sorted.index(current_snarl_id) - 1] if (snarl_ids_sorted.index(current_snarl_id) - 1 >= 0) else -1 
                    left_snarl_end_node_id = -1
                    if (
                            self.snarl_to_anchors_dictionary.get(left_snarl_id) != None
                            and len(self.snarl_to_anchors_dictionary[left_snarl_id]) > 0
                        ):
                        # Get left snarl's end node. We need this to later check if our current snarl could be extended in the left direction 
                        # or not, i.e. if the left snarl has already extended it's end boundary, we need that information.
                        left_snarl_end_node_id =  max(self.snarl_to_anchors_dictionary[left_snarl_id][0][0].id, self.snarl_to_anchors_dictionary[left_snarl_id][0][-1].id)
                        print(f"Snarl on the left is {left_snarl_id}, and it's end node is {left_snarl_end_node_id}")
                    
                    left_degree = self.graph.get_degree(current_snarl_start_handle, go_left_bool)
                    print(f"..left_degree of {current_snarl_start_id} is {left_degree}")
                    if (
                        left_degree == 2
                        and left_snarl_end_node_id != -1
                        and left_snarl_end_node_id == current_snarl_start_id
                        and self.snarl_to_anchors_dictionary[left_snarl_id][0].bp_occupied_end_node + self.snarl_to_anchors_dictionary[current_snarl_id][0].bp_occupied_start_node == self.graph.get_length(current_snarl_start_handle)
                    ):
                        # try merging to one direction
                        print(f"    ..Trying _extending_anchors_by_merging in left")
                        current_snarl_anchors, snarl_ids_list_idx = self._extending_anchors_by_merging(snarl_ids_sorted, snarl_ids_list_idx, current_snarl_id, left_snarl_id, current_snarl_anchors, extend_left=extend_left, anchors_to_discard=anchors_to_remove, snarl_orientation=snarl_orientation, merging_round=merging_round)
                        print(f"    ..#anchors returned after merging snarls {current_snarl_id} and {left_snarl_id}: ", len(current_snarl_anchors))
                        print(f"    ..new snarl id after merging is: {current_snarl_anchors[0].snarl_id}")
                        if current_snarl_anchors[0].snarl_id != current_snarl_id:
                            valid_anchors.extend([[anchor_i, anchor_i.bp_matched_reads] for anchor_i in current_snarl_anchors])
                    min_anchor_length_in_snarl = min([anchor.basepairlength for anchor in current_snarl_anchors])
                    if min_anchor_length_in_snarl >= MIN_ANCHOR_LENGTH:
                        break

                    extend_left = not extend_left
                    current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]
                    current_snarl_end_id = max(current_snarl_anchors[0][0].id, current_snarl_anchors[0][-1].id)
                    print(f"Current snarl end node is {current_snarl_end_id}")
                    current_snarl_end_handle = self.graph.get_handle(current_snarl_end_id)
                    
                    right_snarl_id = snarl_ids_sorted[snarl_ids_sorted.index(current_snarl_id) + 1] if (snarl_ids_sorted.index(current_snarl_id) + 1 < len(snarl_ids_sorted)) else -1 
                    right_snarl_start_node_id = 100000000000000
                    if (
                        self.snarl_to_anchors_dictionary.get(right_snarl_id) != None
                        and len(self.snarl_to_anchors_dictionary[right_snarl_id]) > 0
                    ):
                        right_snarl_start_node_id = min(self.snarl_to_anchors_dictionary[right_snarl_id][0][0].id, self.snarl_to_anchors_dictionary[right_snarl_id][0][-1].id) 
                        print(f"Snarl on the right is {right_snarl_id}, and it's start node is {right_snarl_start_node_id}")

                    right_degree = self.graph.get_degree(current_snarl_end_handle, not go_left_bool)
                    print(f"..right_degree of {current_snarl_end_id} is {right_degree}")
                    if (
                        right_degree == 2
                        and right_snarl_start_node_id != 100000000000000
                        and right_snarl_start_node_id == current_snarl_end_id
                        and self.snarl_to_anchors_dictionary[right_snarl_id][0].bp_occupied_start_node + self.snarl_to_anchors_dictionary[current_snarl_id][0].bp_occupied_end_node == self.graph.get_length(current_snarl_end_handle)
                    ):
                        # current_snarl_id is fetched from list again, as it might have been updated in left-extension
                        current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]
                        print(f"    ..Trying _extending_anchors_by_merging in right")
                        current_snarl_anchors, snarl_ids_list_idx = self._extending_anchors_by_merging(snarl_ids_sorted, snarl_ids_list_idx, current_snarl_id, right_snarl_id, current_snarl_anchors, extend_left=extend_left, anchors_to_discard=anchors_to_remove, snarl_orientation=snarl_orientation, merging_round=merging_round)
                        print(f"    ..#anchors returned after merging snarls {current_snarl_id} and {right_snarl_id}: ", len(current_snarl_anchors))
                        print(f"    ..new snarl id after merging is: {current_snarl_anchors[0].snarl_id}")

                        if current_snarl_anchors[0].snarl_id != current_snarl_id:
                            valid_anchors.extend([[anchor_i, anchor_i.bp_matched_reads] for anchor_i in current_snarl_anchors])
                        print(f"finished extending snarl {current_snarl_id}")
                    min_anchor_length_in_snarl = min([anchor.basepairlength for anchor in current_snarl_anchors])
                current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]

            snarl_ids_list_idx += 1

        # now loop over valid_anchors dict to drop all anchors in anchors_to_remove
        for anchor, reads in valid_anchors:
            if isinstance(anchor.snarl_id, str) and "-" in anchor.snarl_id:
                print(f"snarl {anchor.snarl_id} before adding")
            if anchor not in anchors_to_remove:
                if isinstance(anchor.snarl_id, str) and "-" in anchor.snarl_id:
                    print(f"snarl {anchor.snarl_id} after adding")
                read_info_for_anchor_to_shasta = [read[:4] for read in anchor.bp_matched_reads]
                valid_anchor_extended.append([anchor, read_info_for_anchor_to_shasta])

        return valid_anchor_extended


    def next_handle_iteratee(self, next_boundary):
        self.next_handle_expand_boundary = next_boundary
        # returning False as there is just 1 node connected when the degree is 1.
        return False
    

    def dump_valid_anchors(self, out_file_path, extended_out_file_path, anchor_read_tracking_file_path) -> list:
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
                if len(reads) > MIN_ANCHOR_READS:
                    anchor = self.sentinel_to_anchor[sentinel][id]
                    snarl_id = anchor.snarl_id
                    self.snarl_to_anchor_reads_dictionary[snarl_id].append(len(reads))
                    self.snarl_to_anchors_dictionary[snarl_id].append(anchor)    # stores snarl to anchors mapping for anchor extension

                    anchor_reads = []
                    for read in reads:
                        read[1] = 0 if read[1] else 1
                        anchor_reads.append(read)
                        # Record final coverage
                        self.anchor_coverage.record_final_coverage(f"{anchor!r}", read[0])

                    anchor = self.sentinel_to_anchor[sentinel][id]
                    valid_anchors.append([f"{anchor!r}", anchor_reads])
                    valid_anchors_to_extend.append([anchor, anchor_reads])
        
        dump_to_jsonl(valid_anchors, out_file_path)    

        # extension 
        # self.valid_anchors_extended = self.merge_anchors(valid_anchors_to_extend)   # make sure that it returns serialized anchor object
        self.valid_anchors_extended = self.extend_and_merge_snarls(valid_anchors_to_extend)   # make sure that it returns serialized anchor object
        dump_to_jsonl([[f"{anchor!r}", reads] for anchor, reads in self.valid_anchors_extended], extended_out_file_path)   # also dumping valid_anchors_extended
        dump_to_jsonl(self.anchor_read_tracking_dict, anchor_read_tracking_file_path)    # currently, read drop during snarl merging is not being tracked

        # Save coverage statistics
        self.anchor_coverage.save_coverage_stats(out_file_path + ".coverage.json")
        return valid_anchors


    def dump_anchor_information(self, out_file_path):
        with open(out_file_path, "w") as f:
            print("SNARL_ID\tSENTINEL\tANCHOR\tIS_HETEROZYGOUS")
            for sentinel, anchors_list in self.sentinel_to_anchor.items():
                for anchor in anchors_list:
                    snarl_id = anchor.snarl_id
                    is_heterozygous = (True if (len(self.snarl_to_anchor_reads_dictionary[snarl_id]) > 1) else False)
                    print(f'{snarl_id}\t{sentinel}\t{anchor!r}\t{is_heterozygous}', file=out_file_path)


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

    def dump_dictionary_with_reads_counts(self,out_file_path: str) -> None:
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
        It processes an alignment list (the result of parsing an alignment line) to find anchors in the read associated with the alignment. It updates the anchors dictionary by adding the read to the anchor tuple.
        It:
        1 - walks on the nodes of the path where the read aligns
        2 - if the node where it is standing is a sentinel, checks if there is an anchor that matches the path around the node.
        3 - if it finds one, it then verifies, using the cs tag if the sequence of the reads aligns perfectly (total match) with the portion of the path that is the anchor.
        4 - if so, it adds the read to the list of reads associated to the anchor and stops looking for anchors, as it has to be unique.
        5 - keeps walking until the aligned path ends.

        The read info assigned to the anchor is designed as follows.
        (ReadId, Strand, Begin, End) where ReadId references a read name in the fasta/fastq file, Strand is 0 if the read is oriented as in the fasta/fastq and 1 for the opposite orientation, and Begin, End are the first base of that read in the anchor and the first base in the read following the anchor (that is, End points to one base past the anchor end). Begin, End are coordinates on the read AFTER reverse complementing if Strand is 1. This means that End > Begin always, and the length of the anchor is End - Begin, a positive number.

        Parameters
        ----------
        alignment_l : list
            the list of variables obtained from processing the alignment

        """
        # walk through the nodes keeping the position in the node list
        # and the total length of the nodes
        walked_length = 0
        read_id = alignment_l[READ_POSITION]
        print(f"Processing read {read_id}.....")

        for position, node_id in enumerate(alignment_l[NODE_POSITION]):

            # Verifying that the nodes coming from the alingment are in the graph I am using
            if not self.graph.has_node(node_id):
                print(f"THE NODE {node_id} PRESENT IN THE ALIGNMENT IS NOT IN THE PACKED GRAPH.")
                exit(1)

            node_handle = self.graph.get_handle(node_id)
            length = self.graph.get_length(node_handle)

            anchors = self.sentinel_to_anchor.get(node_id)
            
            if anchors:
                for index, anchor in enumerate(anchors):
                    print(f"For read {read_id}, checking path concordance for anchor {anchor!r}")

                    # an anchor is a list tuple of a list of node handles
                    # and a counter set to 0 at the beginning

                    # scan backward to check that the alignment corresponds to the anchor.
                    alignment_matches_anchor, walk_start, walk_end, relative_strand, walk_start_for_cs_matching, walk_end_for_cs_matching = (
                        verify_path_concordance(
                            position,
                            node_id,
                            alignment_l[NODE_POSITION],
                            alignment_l[ORIENTATION_POSITION],
                            anchor,
                            walked_length
                        )
                    )
                    # if len(anchor) == 2:
                    #     bp_passed -= length
                    #     bp_to_pass += length
                    if alignment_matches_anchor:
                        print(f" {anchor!r} path matched")
                        # if alignment_l[READ_P] == "m64012_190920_173625/50988488/ccs":
                        self.reads_matching_anchor_path += 1
                        # anchor.path_matched_reads.append(read_id)
                        x = (
                            walk_start,
                            walk_end,
                            alignment_l[CIGAR_POSITION],
                            alignment_l[START_POSITION],
                            alignment_l[END_POSITION],
                            walk_start_for_cs_matching,
                            walk_end_for_cs_matching
                        )

                        is_aligning, read_start, read_end, match_limit, cs_start_pos, cs_end_pos = (
                            verify_sequence_agreement(*x)
                        )
                        # If paths is correct:
                        # I need to append the read info to the anchor.
                        # I need read start and read end of the anchor and the orientation of the read
                        if (debug_file):
                            print(f"{read_id},{repr(anchor)},{alignment_matches_anchor},{is_aligning},{match_limit},{cs_start_pos},{cs_end_pos}", file=debug_file)
                        if is_aligning:
                            print(f" {anchor!r} bp matched")
                            # if alignment_l[READ_P] == "m64012_190920_173625/50988488/ccs":
                            self.reads_matching_anchor_sequence += 1                            
                            # if not (alignment_l[STRAND_POSITION]):
                            # TODO: Better relative strand calculation. For first read in anchor, store 0 strand and coordinates. Compute the alignment string.
                            # For next read, if the string is same as previous, then strand = 0, else check if it's reverse complement, then strand = 1. If nothing, then report.
                            if not (relative_strand):
                                tmp = read_start
                                read_start = alignment_l[R_LEN_POSITION] - read_end
                                read_end = alignment_l[R_LEN_POSITION] - tmp

                            # strand = 0 if alignment_l[STRAND_POSITION] else 1
                            strand = 0 if relative_strand else 1
                            anchor.bp_matched_reads.append([alignment_l[READ_POSITION], strand, read_start, read_end, match_limit, cs_start_pos, cs_end_pos])
                            self.anchor_reads_dict[node_id][index].append(
                                [
                                    alignment_l[READ_POSITION],
                                    # alignment_l[STRAND_POSITION],
                                    relative_strand,
                                    read_start,
                                    read_end
                                ]
                            )
                            # Record initial coverage
                            self.anchor_coverage.record_initial_coverage(f"{anchor!r}", alignment_l[READ_POSITION])
                            self.sentinel_to_anchor[node_id][index].add_sequence()
                            # found, no need to check in other anchors
                            break
            # adding to the walked length the one of the node I just passed

            walked_length += length


    
    # def extend_valid_anchors(self, valid_anchors):
    #     # loop over anchors, check if extension required
    #     for anchor, reads in valid_anchors:
    #         anchor_orientation = 
    #         if anchor.baseparilength < MIN_ANCHOR_LENGTH:
    #             start_node, end_node = anchor[0], anchor[-1]
    #             go_left = True
    #             can_go_left, can_go_right = True, True
    #             curr_length = anchor.baseparilength
    #             while (curr_length < MIN_ANCHOR_LENGTH) and (can_go_left or can_go_right):
    #                 # TODO: Update go_left based on loss function 
    #                 extending_node = anchor[0]
    #                 if not go_left:
    #                     extending_node = anchor[-1]
    #                 current_handle = self.graph.get_handle(extending_node.id)
    #                 if extending_node.bp_used_in_extension != 0:
    #                     if go_left:
    #                         can_go_left = False
    #                     else:
    #                         can_go_right = False
    #                     continue

    #                 degree = self.graph.get_degree(current_handle, go_left)
    #                 if degree == 1:
    #                     self.graph.follow_edges(current_handle, go_left, self.next_handle_iteratee)
    #                     if self.next_handle_expand_boundary is None or self.next_handle_expand_boundary == current_handle:
    #                         if go_left:
    #                             can_go_left = False
    #                         else:
    #                             can_go_right = False
    #                         continue
                        
    #                     # now extend 1 node
    #                     if go_left:
    #                         anchor.insert()
    #                     else:
    #                         anchor.
    #                 go_left = not go_left
    #                 pass
        

    #         if degree == 1:
    #             print(f"inside extension, current 1-degree node being checked: {current_handle}")
    #             self.graph.follow_edges(current_handle, go_left_bool, self.next_handle_iteratee)
    #             if self.next_handle_expand_boundary is None or self.next_handle_expand_boundary == current_handle:
    #                 break
    #             else:
    #                 self.anchor_length_occupied += self.graph.get_length(current_handle) - (self.graph.get_length(current_handle) // 2) + (self.graph.get_length(self.next_handle_expand_boundary) // 2)
    #                 nodes_inside_snarl.append(self.graph.get_id(current_handle))
    #                 # nodes_inside_snarl.append(self.graph.get_id(self.next_handle_expand_boundary))
    #                 current_handle = self.next_handle_expand_boundary

    #         else:
    #             break



    #     return valid_anchors




# class AlignSnarlBoundary:

#     def __init__(self) -> None:
#         # useful initialization objects
#         self.graph = PackedGraph()
#         self.boundary_: dict = dict()
#         self.anchor_reads_dict: dict = dict()
#         self.reads_matching_anchor_path: int = 0
#         self.reads_matching_anchor_sequence: int = 0

#     def build(self, dict_path: str, packed_graph_path: str) -> None:

#         # loading dictionary
#         with open(dict_path, 'rb') as in_f:
#             self.sentinel_to_anchor = pickle.load(in_f)

#         #loading packedgraph
#         self.graph.deserialize(packed_graph_path)

#         # initializing output dictionary
#         for sentinel, anchors in self.sentinel_to_anchor.items():
#             self.anchor_reads_dict[sentinel] = [[] for _ in range(len(anchors))]


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

    # INITIALIZING A LIST OF SENTINEL NODE LENGTHS TO 0.
    sentinel_list = anchor.get_sentinels()
    sentinel_bp_consumed_list = [0] * len(sentinel_list)

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
    node_orientations_in_anchor = []

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
        node_orientations_in_anchor.append(alignment_orientation_list[alignment_pos])
        
        #ADDING THE BASEPAIR LENGTHS
        basepairs_consumed_list[anchor_pos] = anchor_concordant[anchor_pos].length

        #ADDING SENTINEL NODE LENGTHS if anchor_pos points to one of the sentinel nodes 
        if anchor_concordant[anchor_pos] in sentinel_list:
            sentinel_bp_consumed_list[sentinel_pos] = anchor_concordant[anchor_pos].length
            sentinel_pos += 1

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
    
    # COMPUTING START AND END OF WALK FOR BASEPAIR SEQUENCE AGREEMENT
    start_walk = walked_length - sum(basepairs_consumed_list[0:sentinel_cut]) + basepairs_consumed_list[0] - (0 if (basepairs_consumed_list[0] == 1) else 1)
    end_walk = walked_length + sum(basepairs_consumed_list[sentinel_cut:]) - basepairs_consumed_list[-1] + (0 if (basepairs_consumed_list[-1] == 1) else 1)
    start_walk_for_cs_matching = start_walk - 1
    end_walk_for_cs_matching = end_walk + 1

    # COMPUTING READ RELATIVE STRAND
    if len(anchor._reads) == 0:
        count_positive_orientation_nodes = node_orientations_in_anchor.count(True)
        relative_strand = True if count_positive_orientation_nodes > (len(node_orientations_in_anchor) / 2) else False
        anchor._reads.append((node_orientations_in_anchor, relative_strand))
    else:
        first_read_orientations, first_read_strand = anchor._reads[0]
        is_orientation_matching = (
            node_orientations_in_anchor == first_read_orientations
        )
        if is_orientation_matching:
            relative_strand = first_read_strand
        else:
            relative_strand = not first_read_strand

    return (True, start_walk, end_walk, relative_strand, start_walk_for_cs_matching, end_walk_for_cs_matching)


def verify_sequence_agreement(
    # self,
    anchor_bp_start: int,
    anchor_bp_end: int,
    cs_walk: list,
    start_in_path: int,
    end_in_path: int,
    walk_start_for_cs_matching: int,
    walk_end_for_cs_matching: int
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
    # If anchor overflows the alingment, it is not valid
    if anchor_bp_end > end_in_path or anchor_bp_start < start_in_path or anchor_bp_end < anchor_bp_start:
        return (False, 0, 0, 0, 0, 0)

    walked_in_the_sequence: int = (
        0  # I need this to keep track of anchor position in the sequence
    )
    walked_in_the_path: int = (
        start_in_path  # I need this to keep track of my walk in the path
    )
    allow_seq_diff: bool = (
        True  # I need this to control no variation between anchor and sequence is present. Starting with True, setting to False when walking on anchor coordinates
    )

    # When walking on alingment. Path length is calculated as 'equal + subst + delition'
    # When walking on alingment. Read length is calculated as 'equal+subst+insertion'
    # For the moment, strand can be assumed as +
    total_matched_bps = 0
    for step in cs_walk:

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
            # If I passed on a equal step, it is ok. I set allow_differences to false and go on. But before I check if I have surpassed the end of the anchor. If yes return true.
            total_matched_bps = step[1]
            if walked_in_the_path >= walk_end_for_cs_matching:
                diff_start = walked_in_the_path - anchor_bp_start
                diff_end = walked_in_the_path - anchor_bp_end
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

        # I passed the end of the scan and there was no difference
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
                total_matched_bps - diff_end
            )

        elif walked_in_the_path > end_in_path:
            return (False, 0, 0, 0, 0, 0)
    return (False, 0, 0, 0, 0, 0)
