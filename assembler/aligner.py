import json
from sys import stderr
import pickle
import os.path
from sys import exit
from collections import defaultdict
import copy

from bdsg.bdsg import PackedGraph
from assembler.anchor import Anchor
from assembler.constants import (
    READ_POSITION,
    R_LEN_POSITION,
    STRAND_POSITION,
    START_POSITION,
    END_POSITION,
    NODE_POSITION,
    ORIENTATION_POSITION,
    CIGAR_POSITION,
    READS_DEPTH_HETEROZYGOUS,
    READS_DEPTH_HOMOZYGOUS,
    SNARL_DEPTH_HETEROZYGOUS,
    MIN_ANCHOR_READS,
    MIN_ANCHOR_LENGTH,
    READ_POSITION,
    HET_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING,
    HOMO_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING,
    MIN_READS_REQUIRED_FOR_MERGING,
)

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


    def _extending_anchors_by_merging(self, snarl_ids_sorted_list_up_to_date, snarl_ids_sorted_list_iterator_idx, current_snarl_id, other_snarl_id, current_snarl_anchors, extend_left, anchors_to_discard, snarl_orientation) -> list:
        
        print("current snarl being extended:- ", current_snarl_id)
        print("surrounding snarls:- ", snarl_ids_sorted_list_up_to_date[snarl_ids_sorted_list_iterator_idx-5:min(len(snarl_ids_sorted_list_up_to_date), snarl_ids_sorted_list_iterator_idx+5)])
        
        # calculating if after merging, enough anchors (>=2) will remain 
        cnt_anchors_with_sufficient_read_overlap = 0
        # other_anchor = self.snarl_to_anchors_dictionary[other_snarl_id][0]
        FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING = HET_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING
        if len(self.snarl_to_anchors_dictionary[other_snarl_id]) == 1:    # for homozygous snarls, use homozygous read retain fraction while merging
            FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING = HOMO_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING

        for other_anchor in self.snarl_to_anchors_dictionary[other_snarl_id]:
            for anchor in current_snarl_anchors:
                common_paths = set(anchor.reference_paths_covered).intersection(set(other_anchor.reference_paths_covered))
                read_ids_current_anchor = [read[READ_POSITION] for read in anchor.bp_matched_reads]
                read_ids_other_anchor = [read[READ_POSITION] for read in other_anchor.bp_matched_reads]
                common_reads_ids = set(read_ids_current_anchor).intersection(set(read_ids_other_anchor)) 
                fraction_reads_retained_in_current = len(common_reads_ids) / len(set(read_ids_current_anchor))
                fraction_reads_retained_in_other = len(common_reads_ids) / len(set(read_ids_other_anchor))
                if len(common_paths) > 0 and (len(common_reads_ids) > MIN_READS_REQUIRED_FOR_MERGING) and (fraction_reads_retained_in_current >= FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING) and (fraction_reads_retained_in_other >= FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING):    # meaning we can create an anchor with this combination
                    cnt_anchors_with_sufficient_read_overlap += 1
                    if cnt_anchors_with_sufficient_read_overlap > 1:
                        break
                    # extendable_anchors_in_current_snarl.append(anchor)
            if cnt_anchors_with_sufficient_read_overlap > 1:
                break

        if cnt_anchors_with_sufficient_read_overlap < 2:
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
                fraction_reads_retained_in_current = len(common_reads_ids) / len(set(read_ids_current_anchor))
                fraction_reads_retained_in_other = len(common_reads_ids) / len(set(read_ids_other_anchor))
                if len(common_paths) > 0 and (len(common_reads_ids) > MIN_READS_REQUIRED_FOR_MERGING) and (fraction_reads_retained_in_current >= FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING) and (fraction_reads_retained_in_other >= FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING):    # meaning we can create an anchor with this combination
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
                    # TODO: Fix the new_snarl_id generation
                    new_snarl_id = (str(other_anchor.snarl_id) + "-" + str(new_anchor.snarl_id)) if extend_left else (str(new_anchor.snarl_id) + "-" + str(other_anchor.snarl_id))
                    new_anchor.add_snarl_id(new_snarl_id)
                    new_anchor.compute_bp_length()
                    new_anchor.reference_paths_covered = set(new_anchor.reference_paths_covered).intersection(set(other_anchor.reference_paths_covered))
                    # common_read_ids = set([read[READ_POSITION] for read in anchor.bp_matched_reads]).intersection(set([read[READ_POSITION] for read in other_anchor.bp_matched_reads]))
                    # find start, end of common reads
                    common_bp_matched_reads = {}
                    for read in new_anchor.bp_matched_reads:
                        if read[0] in common_reads_ids:
                            common_bp_matched_reads[read[0]] = read
                    for read in other_anchor.bp_matched_reads:
                        if read[0] in common_reads_ids:
                            unpacked_read_id, unpacked_strand, unpacked_start, unpacked_end = common_bp_matched_reads[read[0]]
                            unpacked_start = min(read[2], unpacked_start)
                            unpacked_end = max(read[3], unpacked_end)
                            # unpacked_strand_orientation = True if unpacked_strand == 0 else 1
                            # if (
                            #     ((extend_left) and (snarl_orientation == unpacked_strand_orientation)) 
                            #     or ((not extend_left) and (snarl_orientation != unpacked_strand_orientation))):
                            #     unpacked_start = read[2]
                            # else:
                            #     unpacked_end = read[3]
                            common_bp_matched_reads[read[0]] = [unpacked_read_id, unpacked_strand, unpacked_start, unpacked_end]
                    
                    new_anchor.bp_matched_reads = list(common_bp_matched_reads.values())    # set(anchor.bp_matched_reads).intersection(set(other_anchor.bp_matched_reads))
                    new_anchors_after_merging.append(new_anchor)
                # else:    # homozygous merging codeblock (old)
                #     # remove anchor
                #     anchors_to_discard.append(anchor)

        # remove all anchors from both current and other snarls, as they are now replaced by new anchors having new snarl name
        anchors_to_discard.extend(current_snarl_anchors)
        anchors_to_discard.extend(self.snarl_to_anchors_dictionary[other_snarl_id])
        # add the new snarl and its anchors in the snarl_to_anchors_dictionary
        new_snarl_id_after_merge = new_anchors_after_merging[0].snarl_id
        self.snarl_to_anchors_dictionary[new_snarl_id_after_merge] = current_snarl_anchors
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


    def merge_anchors(self, valid_anchors: list) -> list:
        """
        * Iterate over shorter anchors, find adjacent snarls (+1/-1). If read drop from one snarl to the other is within the defined threshold,
        then merge the snarls. Get all combinations of anchors (required it belongs to atleast one path) and re-define this as a new anchor,
        with common reads and newly computed bplength, snarl_id and pathnames.

        * If both left, right anchors are available for merging, choose the direction where read coverage drop is minimum
            # get snarl direction (0 -> forward, 1 -> backward)
            # when merging anchors, if orientation of other anchor is opposite to that of current, do the following:-
            # 1.) flip other anchor
            # 2.) if other anchor can be found by extending in same direction as current anchor, then add other anchor to the end of current anchor. else, add to beginning.
            # TODO: for calculating new start and end posn of reads

        """
        valid_anchor_extended = []
        # newsnarls_to_anchors_dictionary = {}
        anchors_to_remove = []   # {(snarl_id, anchor)}

        snarl_orientation = True    # True-> snarl is going in forward dirn, False otherwise.
        snarl_id_1, snarl_id_2 = -1, -1
        for snarl_id, snarl_anchors in self.snarl_to_anchors_dictionary.items():
            if len(snarl_anchors) > 0:
                if snarl_id_1 == -1:
                    snarl_id_1 = snarl_id
                else:
                    snarl_id_2 = snarl_id
                    break
        if snarl_id_1 == -1 or snarl_id_2 == -1:
            return valid_anchors
        # an edge case where this calculation might go wrong is if we have two adjacent snarls having anchors with only 2 nodes 
        if (snarl_id_1 - snarl_id_2) * (self.snarl_to_anchors_dictionary[snarl_id_1][0][1].id - self.snarl_to_anchors_dictionary[snarl_id_2][0][0].id) < 0:    # means opposite direction
            snarl_orientation = False

        snarl_ids_sorted = sorted(list(self.snarl_to_anchors_dictionary.keys()))
        snarl_ids_list_idx = 0    # we need to use this index counter (and can't simply use an iterator) because we will be inserting/deleting the snarl_ids_sorted list on the go
        # Note: Remember to update the snarl_ids_list_idx appropriately when merging snarls
        while snarl_ids_list_idx < len(snarl_ids_sorted):
            current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]
            last_snarl_id = -1
            current_snarl_anchors = self.snarl_to_anchors_dictionary[current_snarl_id]
            min_anchor_length_in_snarl = min([anchor.baseparilength for anchor in current_snarl_anchors])
            # if len(current_snarl_anchors) > 1 and min_anchor_length_in_snarl < MIN_ANCHOR_LENGTH:
            while (
                len(current_snarl_anchors) > 1 
                and (min_anchor_length_in_snarl < MIN_ANCHOR_LENGTH 
                and last_snarl_id != current_snarl_id)
            ):
                current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]
                last_snarl_id = current_snarl_id
                left_snarl_id = snarl_ids_sorted[snarl_ids_sorted.index(current_snarl_id) - 1] if (snarl_ids_sorted.index(current_snarl_id) - 1 > 0) else -1 
                right_snarl_id = snarl_ids_sorted[snarl_ids_sorted.index(current_snarl_id) + 1] if (snarl_ids_sorted.index(current_snarl_id) + 1 < len(snarl_ids_sorted)) else -1 
                extend_left = True
                if (
                    left_snarl_id == -1 
                    or (isinstance(left_snarl_id, str) and ('-' in left_snarl_id))
                ):
                    extend_left = not extend_left
                    tmp = left_snarl_id
                    left_snarl_id = right_snarl_id
                    right_snarl_id = tmp

                if (
                    self.snarl_to_anchors_dictionary.get(left_snarl_id) != None
                    and len(self.snarl_to_anchors_dictionary[left_snarl_id]) > 0
                ):
                    # try extending to one direction
                    current_snarl_anchors, snarl_ids_list_idx = self._extending_anchors_by_merging(snarl_ids_sorted, snarl_ids_list_idx, current_snarl_id, left_snarl_id, current_snarl_anchors, extend_left=extend_left, anchors_to_discard=anchors_to_remove, snarl_orientation=snarl_orientation)
                    min_anchor_length_in_snarl = min([anchor.baseparilength for anchor in current_snarl_anchors])
                    print(f"#anchors returned after merging snarls {current_snarl_id} and {left_snarl_id}: ", len(current_snarl_anchors))
                    print(f"snarl_id: {current_snarl_anchors[0].snarl_id}")
                    if current_snarl_anchors[0].snarl_id != current_snarl_id:
                        valid_anchors.extend([[anchor_i, anchor_i.bp_matched_reads] for anchor_i in current_snarl_anchors])
                if min_anchor_length_in_snarl >= MIN_ANCHOR_LENGTH:
                    break
                extend_left = not extend_left
                # trying extending to the other direction, only if enough anchors left in current snarl after first merging
                if (
                    self.snarl_to_anchors_dictionary.get(right_snarl_id) != None
                    and len(self.snarl_to_anchors_dictionary[right_snarl_id]) > 0
                ):
                    # current_snarl_id is fetched from list again, as it might have been updated in left-extension
                    current_snarl_id = snarl_ids_sorted[snarl_ids_list_idx]
                    current_snarl_anchors, snarl_ids_list_idx = self._extending_anchors_by_merging(snarl_ids_sorted, snarl_ids_list_idx, current_snarl_id, right_snarl_id, current_snarl_anchors, extend_left=extend_left, anchors_to_discard=anchors_to_remove, snarl_orientation=snarl_orientation)
                    min_anchor_length_in_snarl = min([anchor.baseparilength for anchor in current_snarl_anchors])
                    print(f"#anchors returned after merging snarls {current_snarl_id} and {right_snarl_id}: ", len(current_snarl_anchors))
                    print(f"snarl_id: {current_snarl_anchors[0].snarl_id}")
                    if current_snarl_anchors[0].snarl_id != current_snarl_id:
                        valid_anchors.extend([[anchor_i, anchor_i.bp_matched_reads] for anchor_i in current_snarl_anchors])
                print(f"finished extending snarl {current_snarl_id}")
            snarl_ids_list_idx += 1
                    

        # now loop over valid_anchors dict to drop all anchors in anchors_to_remove
        for anchor, reads in valid_anchors:
            if isinstance(anchor.snarl_id, str) and "-" in anchor.snarl_id:
                print(f"snarl {anchor.snarl_id} before adding")
            if anchor not in anchors_to_remove:
                if isinstance(anchor.snarl_id, str) and "-" in anchor.snarl_id:
                    print(f"snarl {anchor.snarl_id} after adding")
                valid_anchor_extended.append([anchor, anchor.bp_matched_reads])

        return valid_anchor_extended


    
    def dump_valid_anchors(self, out_file_path, extended_out_file_path) -> list:
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
            for id, reads in enumerate(self.anchor_reads_dict[sentinel]):
                if len(reads) > MIN_ANCHOR_READS:
                    anchor = self.sentinel_to_anchor[sentinel][id]
                    snarl_id = anchor.snarl_id
                    self.snarl_to_anchor_reads_dictionary[snarl_id].append(len(reads))
                    self.snarl_to_anchors_dictionary[snarl_id].append(anchor)    # stores snarl to anchors mapping for anchor extension

        for sentinel in self.anchor_reads_dict:
            for id, reads in enumerate(self.anchor_reads_dict[sentinel]):
                sentinel_anchor = []
                snarl_id = self.sentinel_to_anchor[sentinel][id].snarl_id

                if (len(self.snarl_to_anchor_reads_dictionary[snarl_id]) > 1 and len(reads) >= READS_DEPTH_HETEROZYGOUS and sum(self.snarl_to_anchor_reads_dictionary[snarl_id]) >= SNARL_DEPTH_HETEROZYGOUS) or (len(reads) >= READS_DEPTH_HOMOZYGOUS):
                    for read in reads:
                        read[1] = 0 if read[1] else 1
                        sentinel_anchor.append(read)
                    if len(sentinel_anchor) > MIN_ANCHOR_READS:
                        anchor = self.sentinel_to_anchor[sentinel][id]
                        valid_anchors.append([f"{anchor!r}", sentinel_anchor])
                        valid_anchors_to_extend.append([anchor, sentinel_anchor])
        
        dump_to_jsonl(valid_anchors, out_file_path)    

        # extension 
        self.valid_anchors_extended = self.merge_anchors(valid_anchors_to_extend)   # make sure that it returns serialized anchor object
        dump_to_jsonl([[f"{anchor!r}", reads] for anchor, reads in self.valid_anchors_extended], extended_out_file_path)   # also dumping valid_anchors_extended

    def dump_anchor_information(self, out_file_path: str):
        with open(out_file_path, "w") as f:
            print("SNARL_ID\tSENTINEL\tANCHOR\tIS_HETEROZYGOUS")
            for sentinel in self.sentinel_to_anchor:
                for anchor in self.sentinel_to_anchor[sentinel]:
                    is_heterozygous = True if len(self.snarl_to_anchor_reads_dictionary[anchor.snarl_id]) > 1 else False
                    print(f'{anchor.snarl_id}\t{sentinel}\t{anchor!r}\t{is_heterozygous}')
        # valid_anchors = []
        # for sentinel in self.anchor_reads_dict:
        #     for id, reads in enumerate(self.anchor_reads_dict[sentinel]):
        #         sentinel_anchor = []
        #         if len(reads) > READS_DEPTH:
        #             for read in reads:
        #                 read[1] = 0 if read[1] else 1
        #                 sentinel_anchor.append(read)
        #         if len(sentinel_anchor) > 0:
        #             anchor = self.sentinel_to_anchor[sentinel][id]
        #             valid_anchors.append([f"{anchor!r}", sentinel_anchor])

        # self.dump_to_jsonl(valid_anchors, out_file_path)

    def print_extended_anchor_info(self, out_f) -> None:
        with open(out_f, "w") as f:
            print(f"Sentinel_node\tsnarl_id\tAnchor_length\tAnchor_pos_in_ref_path\tAnchor_path\tAnchor_nodes_copypaste_bandage\tPaths_associated_with_anchor\tbp_matched_reads",file=f)
            for anchor, _ in self.valid_anchors_extended:
                print(
                    f"{anchor.get_sentinel_id()}\t{anchor.snarl_id}\t{anchor.baseparilength}\t{anchor.genomic_position}\t{anchor!r}\t{anchor.bandage_representation()}\t{anchor.get_reference_paths()}\t{len([x[0] for x in anchor.bp_matched_reads])}",
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
                    print(f"For read {read_id}, checking path and sequence concordance for anchor {anchor!r}")

                    # an anchor is a list tuple of a list of node handles
                    # and a counter set to 0 at the beginning

                    # scan backward to check that the alignment corresponds to the anchor.
                    alignment_matches_anchor, walk_start, walk_end, relative_strand = (
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
                            alignment_l[END_POSITION]
                        )

                        is_aligning, read_start, read_end = (
                            verify_sequence_agreement(*x)
                        )
                        # If paths is correct:
                        # I need to append the read info to the anchor.
                        # I need read start and read end of the anchor and the orientation of the read
                        if (debug_file):
                            print(f"{read_id},{repr(anchor)},{alignment_matches_anchor},{is_aligning}", file=debug_file)
                        if is_aligning:
                            print(f" {anchor!r} bp matched")
                            # if alignment_l[READ_P] == "m64012_190920_173625/50988488/ccs":
                            self.reads_matching_anchor_sequence += 1                            
                            # if not (alignment_l[STRAND_POSITION]):
                            if not (relative_strand):
                                tmp = read_start
                                read_start = alignment_l[R_LEN_POSITION] - read_end
                                read_end = alignment_l[R_LEN_POSITION] - tmp

                            # strand = 0 if alignment_l[STRAND_POSITION] else 1
                            strand = 0 if relative_strand else 1
                            anchor.bp_matched_reads.append([alignment_l[READ_POSITION], strand, read_start, read_end])
                            self.anchor_reads_dict[node_id][index].append(
                                [
                                    alignment_l[READ_POSITION],
                                    # alignment_l[STRAND_POSITION],
                                    relative_strand,
                                    read_start,
                                    read_end
                                ]
                            )
                            self.sentinel_to_anchor[node_id][index].add_sequence()
                            # found, no need to check in other anchors
                            break
            # adding to the walked length the one of the node I just passed

            walked_length += length

    # def next_handle_iteratee(self, next_boundary):
    #     self.next_handle_expand_boundary = next_boundary
    #     # returning False as there is just 1 node connected when the degree is 1.
    #     return False
    
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
        json.dump(object, f, ensure_ascii=False)

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
    # The “cut” value tells the function how far from the start of the anchor the sentinel is located.
    sentinel_cut = (
        (len(anchor) - 1 - sentinel_position)
        if not concordance_orientation
        else sentinel_position
    )

    # POSITION OF THE ALIGNMENT AT THE BEGINNING OF THE ANCHOR. IF < 0 OR GREATER THAN ALIGNMENT NODES, EXIT.
    alignment_pos = alignment_position - sentinel_cut
    if alignment_pos < 0 or alignment_pos >= len(alignment_node_id_list):
        return (False, 0, 0, -1)
    
    # INITIALZING A LIST WITH ANCHOR LENGTH TO ZERO. TO KEEP TRACK OF THE BASEPAIRS CONSUMED
    basepairs_consumed_list = [0] * len(anchor)

    # TO SIMPLIFY OPERATIONS, IF THE ANCHOR IS REVERSED COMPARED TO THE PATH, REVERT THE ANCHOR SO SCANNING IS EASIER
    anchor_concordant = anchor[::-1] if not concordance_orientation else anchor[:]
    
    # POSITION IN SCANNING THE ANCHOR
    anchor_pos = 0

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
            return (False, 0, 0, -1)

        # Store read orientations w.r.t anchor nodes in path
        node_orientations_in_anchor.append(alignment_orientation_list[alignment_pos])
        
        #ADDING THE BASEPAIR LENGTHS
        basepairs_consumed_list[anchor_pos] = anchor_concordant[anchor_pos].length

        # INCREASING POSITION COUNTER
        anchor_pos += 1
        alignment_pos += 1

    if anchor_pos < len(anchor_concordant):
        # didn't finish walking the entire anchor, probably because of alignment_pos < len(alignment_node_id_list)
        return (
            False,
            0,
            0,
            -1
        )
    
    # COMPUTING START AND END OF WALK FOR BASEPAIR SEQUENCE AGREEMENT (COVERTED TO CEIL DIV)
    start_walk = walked_length - sum(basepairs_consumed_list[0:sentinel_cut]) + (basepairs_consumed_list[0] + 1) // 2
    end_walk = walked_length + sum(basepairs_consumed_list[sentinel_cut:]) - (basepairs_consumed_list[-1] + 1) // 2

    # COMPUTING READ RELATIVE STRAND
    count_positive_orientation_nodes = node_orientations_in_anchor.count(True)
    relative_strand = True if count_positive_orientation_nodes > (len(node_orientations_in_anchor) / 2) else False

    return (True, start_walk, end_walk, relative_strand)


def verify_sequence_agreement(
    # self,
    anchor_bp_start: int,
    anchor_bp_end: int,
    cs_walk: list,
    start_in_path: int,
    end_in_path: int,
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
        return (False, 0, 0)

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
                return (False, 0, 0)
            # If I passed on a equal step, it is ok. I set allow_differences to false and go on. But before I check if I have surpassed the end of the anchor. If yes return true.
            if walked_in_the_path >= anchor_bp_end:
                diff_start = walked_in_the_path - anchor_bp_start
                diff_end = walked_in_the_path - anchor_bp_end
                # I add a + 1 in the read_end position because of Shasta requirement that the interval is open at the end. The end id in the sequence is of the first nucleotide after the anchor
                # TODO: Currently end_node_pos causes a gap when anchor end node is even #base-pairs.
                return (
                    True,
                    walked_in_the_sequence - diff_start,
                    walked_in_the_sequence - diff_end,
                )
            else:
                allow_seq_diff = False  # go to the next step

        # Walking in the anchor section and found a diff
        elif not (allow_seq_diff) and step[0] != ":":
            return (False, 0, 0)

        # I passed the end of the scan and there was no difference
        elif walked_in_the_path >= anchor_bp_end:
            diff_start = walked_in_the_path - anchor_bp_start
            diff_end = walked_in_the_path - anchor_bp_end
            # I add a + 1 in the read_end position because of Shasta requirement that the interval is open at the end. The end id in the sequence is of the first nucleotide after the anchor
            return (
                True,
                walked_in_the_sequence - diff_start,
                walked_in_the_sequence - diff_end,
            )

        elif walked_in_the_path > end_in_path:
            return (False, 0, 0)
    return (False, 0, 0)

    