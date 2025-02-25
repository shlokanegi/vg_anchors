import json
from sys import stderr
import pickle
import os.path
from sys import exit

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
    READS_DEPTH,
)

# read_compare = 'm64015_190920_185703/40699337/ccs' #'m64012_190921_234837/35261139/ccs'#'m64012_190921_234837/96012404/ccs'

class AlignAnchor:

    def __init__(self) -> None:
        # useful initialization objects
        self.graph = PackedGraph()
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
    
    def dump_valid_anchors(self, out_file_path) -> list:
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
        for sentinel in self.anchor_reads_dict:
            for id, reads in enumerate(self.anchor_reads_dict[sentinel]):
                sentinel_anchor = []
                if len(reads) > READS_DEPTH:
                    for read in reads:
                        read[1] = 0 if read[1] else 1
                        sentinel_anchor.append(read)
                if len(sentinel_anchor) > 0:
                    anchor = self.sentinel_to_anchor[sentinel][id]
                    valid_anchors.append([f"{anchor!r}", sentinel_anchor])

        self.dump_to_jsonl(valid_anchors, out_file_path)

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


    def dump_to_jsonl(self, object, out_file_path: str):
        """
        It dumps the object to json structure.

        Parameters
        ----------
        valid_anchors : list
            the list containing lists of anchors for each sentinel.
        """
        with open(out_file_path, "w", encoding="utf-8") as f:
            json.dump(object, f, ensure_ascii=False)



    def processGafLine(self, alignment_l: list, debug_file):
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

                    # an anchor is a list tuple of a list of node handles
                    # and a counter set to 0 at the beginning

                    # scan backward to check that the alignment corresponds to the anchor.
                    alignment_matches_anchor, walk_start, walk_end = (
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
                        # if alignment_l[READ_P] == "m64012_190920_173625/50988488/ccs":
                        self.reads_matching_anchor_path += 1
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
                        print(f"{read_id},{repr(anchor)},{alignment_matches_anchor},{is_aligning}", file=debug_file)
                        if is_aligning:
                            # if alignment_l[READ_P] == "m64012_190920_173625/50988488/ccs":
                            self.reads_matching_anchor_sequence += 1 
                            if not (alignment_l[STRAND_POSITION]):
                                tmp = read_start
                                read_start = alignment_l[R_LEN_POSITION] - read_end
                                read_end = alignment_l[R_LEN_POSITION] - tmp
                            
                            self.anchor_reads_dict[node_id][index].append(
                                [
                                    alignment_l[READ_POSITION],
                                    alignment_l[STRAND_POSITION],
                                    read_start,
                                    read_end,
                                ]
                            )
                            self.sentinel_to_anchor[node_id][index].add_sequence()
                            # found, no need to check in other anchors
                            break
            # adding to the walked length the one of the node I just passed

            walked_length += length

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
        return (False, 0, 0)
    
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
            return (False, 0, 0)
        
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
        )
    
    # COMPUTING START AND END OF WALK FOR BASEPAIR SEQUENCE AGREEMENT
    start_walk = walked_length - sum(basepairs_consumed_list[0:sentinel_cut]) + (basepairs_consumed_list[0] + 1) // 2
    end_walk = walked_length + sum(basepairs_consumed_list[sentinel_cut:]) - (basepairs_consumed_list[-1] + 1) // 2

    return (True, start_walk, end_walk)

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
        True  # I need this to control no variation beteen anchor and sequence is present. Starting with True, setting to False when walking on anchor coordinates
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

    
