import json
from sys import stderr
import pickle
import os.path

from bdsg.bdsg import PackedGraph
from assembler.constants import (
    READ_P,
    R_LEN_P,
    STRAND_P,
    START_P,
    END_P,
    NODE_P,
    ORIENT_P,
    CS_P,
    READS_DEPTH,
)


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
            self.anchor_reads_dict[sentinel] = [[] for x in range(len(anchors))]

        # for sentinel in self.sentinel_to_anchor:
        #     print(f"S_T_A {sentinel} = {self.sentinel_to_anchor[sentinel]}")
        #     break


    def processGafLine(self, alignment_l: list):
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

        for position, node_id in enumerate(alignment_l[NODE_P]):
            anchors = self.sentinel_to_anchor.get(node_id)
            # print(f"anchors: {anchors}", flush=True)
            if anchors:
                for index, anchor in enumerate(anchors):

                    # an anchor is a list tuple of a list of node handles
                    # and a counter set to 0 at the beginning
                    # anchor = anchor_capsule[0]
                    # locate the position of the sentinel
                    #print(anchor, flush=True)
                    sentinel_p = next(
                        p
                        for p, a in enumerate(anchor)
                        if a.id == node_id
                    )
                    # print(f"Sentinel {node_id} at position {sentinel_p} out of {len(anchor)}")
                    sentinel_orientation = (
                        True if anchor[sentinel_p].orientation else False
                    )
                    concordance_orientation = (
                        sentinel_orientation == alignment_l[ORIENT_P][position]
                    )

                    # scan backward to check that the alignment corresponds to the anchor.
                    alignment_matches_anchor, bp_passed, bp_to_pass = (
                        self.verify_path_concordance(
                            position,
                            sentinel_p,
                            concordance_orientation,
                            alignment_l[NODE_P],
                            alignment_l[ORIENT_P],
                            anchor,
                        )
                    )

                    if alignment_matches_anchor:
                        # if alignment_l[READ_P] == "m64012_190920_173625/50988488/ccs":
                        self.reads_matching_anchor_path += 1
                        x = (
                            walked_length - bp_passed,
                            walked_length + bp_to_pass,
                            alignment_l[CS_P],
                            alignment_l[STRAND_P],
                            alignment_l[START_P],
                            alignment_l[END_P],
                            alignment_l[READ_P],
                        )

                        is_aligning, read_start, read_end = (
                            self.verify_sequence_agreement(*x)
                        )
                        # If paths is correct:
                        # I need to append the read info to the anchor.
                        # I need read start and read end of the anchor and the orientation of the read
                        if is_aligning:
                            # if alignment_l[READ_P] == "m64012_190920_173625/50988488/ccs":
                            self.reads_matching_anchor_sequence += 1 
                            if not (alignment_l[ORIENT_P][position]):
                                tmp = read_start
                                read_start = alignment_l[R_LEN_P] - read_end
                                read_end = alignment_l[R_LEN_P] - tmp
                            self.anchor_reads_dict[node_id][index].append(
                                [
                                    alignment_l[READ_P],
                                    alignment_l[ORIENT_P][position],
                                    read_start,
                                    read_end,
                                ]
                            )
                            # found, no need to check in other anchors
                            break
                        #else:
                            # if alignment_l[READ_P] == "m64012_190920_173625/50988488/ccs":
                            #     print(f"Not aligning read {alignment_l[READ_P]}. Ranges are {walked_length - bp_passed},{walked_length + bp_to_pass}", file=stderr)
            # adding to the walked length the one of the node I just passed
            node_handle = self.graph.get_handle(node_id)
            length = self.graph.get_length(node_handle)

            walked_length += length

    def verify_path_concordance(
        self,
        alignment_position: int,
        sentinel_position: int,
        concordance_orientation: bool,
        alignment_node_id_list: list,
        alignment_orientation_list: list,
        anchor: list
    ) -> list:
        """
        It verifies that the path around the node where the process_alignment function is standing matches the anchor.
        If so, it returns True and returns how many base pairs before and after the start of the sentinel node the sequence alignment ha to be perfect match to validate the anchor.

        Parameters
        ----------
        alignment_position: int
            The position of the current sentinel node being evaluated by the process_alignment function. Is serves as beginning position to locate the sentinel node in the alignment node list.
        sentinel_position: int
            The position in the anchor of the sentinel node.
        concordance_orientation: bool
            True if the orientation of the sentinel node is the same between anchor and aligned path, else False
        alignment_node_id_list: list
            The list of node_id corresponding to the alignment path
        alignment_orientation_list: list
            The list of node orientations correspoding to the nodes in alignment_node_id_list
        anchor: list
            The anchor list

        Returns
        -------
        bool
            True if iteration has to continue else False
        already_walked: int
            The basepairs between the start of the sentinel node and the start of the anchor
        to_walk: int
            The basepairs between the start of the sentinel node and the end of the anchor

        """

        # if they concorde
        bp_to_walk = [0] * len(anchor)
        sentinel_cut = (
            (len(anchor) - 1 - sentinel_position)
            if not concordance_orientation
            else sentinel_position
        )

        anchor_c = anchor[::-1] if not concordance_orientation else anchor[:]
        anchor_pos = 0
        alignment_pos = alignment_position - sentinel_cut

        if alignment_pos < 0 or alignment_pos >= len(alignment_node_id_list):
            return (False, 0, 0)

        while anchor_pos < len(anchor_c) and alignment_pos < len(
            alignment_node_id_list
        ):
            # check node_id is the same
            if (
                alignment_node_id_list[alignment_pos]
                != anchor_c[anchor_pos].id
            ) or (
                concordance_orientation
                != (
                    alignment_orientation_list[alignment_pos]
                    == anchor_c[anchor_pos].orientation
                )
            ):
                return (False, 0, 0)
            bp_to_walk[anchor_pos] = anchor_c[anchor_pos].length
            anchor_pos += 1
            alignment_pos += 1

        if anchor_pos < len(anchor_c):
            # didn't finish walking the entire anchor. Should be not necessary
            return (
                False,
                0,
                0,
            )

        if len(anchor_c) > 2:
            already_walked = sum(bp_to_walk[1:sentinel_cut]) + bp_to_walk[0] // 2
            to_walk = sum(bp_to_walk[sentinel_cut:-1]) + bp_to_walk[-1] // 2
        else:
            already_walked = bp_to_walk[0] // 2
            to_walk = bp_to_walk[-1] // 2

        return (True, already_walked, to_walk)

    def verify_sequence_agreement(
        self,
        anchor_bp_start: int,
        anchor_bp_end: int,
        cs_walk: list,
        seq_strand: bool,
        start_in_path: int,
        end_in_path: int,
        read_id: str
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
        seq_strand: bool
            Not used for now as all alignments are +? Probably due to the sample I am using or some implementation?
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
        # print(seq_strand)
        # If anchor overflows the alingment, it is not valid
        if anchor_bp_end > end_in_path or anchor_bp_start < start_in_path:
            # if read_id == "m64012_190920_173625/50988488/ccs":
            #     print(f"Mismatch because anchor_bp_end {anchor_bp_end} > end_in_path {end_in_path} or anchor_bp_start {anchor_bp_start} < start_in_path {start_in_path}",file=stderr)
            
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
                    # if read_id == "m64012_190920_173625/50988488/ccs":
                    #     print(f"step {step[0]} is mismatch. Walked in path: {walked_in_the_path}, anchor_start: {anchor_bp_start}.",file=stderr)
                    return (False, 0, 0)
                # If I passed on a equal step, it is ok. I set allo_differences to false and go on. But before I check if I have surpassed the end of the anchor. If yes return true.
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
                # if read_id == "m64012_190920_173625/50988488/ccs":
                #         print(f"step {step[0]} is mismatch. Walked in path: {walked_in_the_path}, anchor_start: {anchor_bp_start}. I am walking in the alingment, end: {anchor_bp_end}",file=stderr)
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
                # if read_id == "m64012_190920_173625/50988488/ccs":
                #         print(f"Surpassed the end of the alignment. Walked in path: {walked_in_the_path}, anchor_start: {anchor_bp_start}. I am walking in the alingment, end: {anchor_bp_end}, end in path {end_in_path}",file=stderr)
                return (False, 0, 0)
        # if read_id == "m64012_190920_173625/50988488/ccs":
        #     print(f"I am at the end and did not take any other return. Walked in path: {walked_in_the_path}, anchor_start: {anchor_bp_start}. I am walking in the alingment, end: {anchor_bp_end}, end in path {end_in_path}",file=stderr)
        return (False, 0, 0)

    def get_anchor_size(self, anchor: list) -> int:
        """
        It returns the size of an anchor. Copy of the one in the anchor dictionaty generation but I need the packedgraph stored in the object to do so.

        Parameters
        ----------
        anchor: list
            the anchor list

        Returns
        -------
        anchor_size: int
            The size in bp of the anchor
        """
        # if a node has odd size, it is goind to be divided as [l,l+1].
        # so for the boundary nodes, in case they are odds, I will take
        # the last l + 1 nucleotides of the starting node
        # the first l nucleotides of the last node
        anchor_size = (
            anchor[0].length // 2
            + anchor[-1].length // 2
        )

        for node_handle in anchor[1:-1]:
            anchor_size += node_handle.length

        return anchor_size

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
            for reads in self.anchor_reads_dict[sentinel]:
                sentinel_anchor = []
                if len(reads) > READS_DEPTH:
                    for read in reads:
                        read[1] = 0 if read[1] else 1
                        sentinel_anchor.append(read)
                if len(sentinel_anchor) > 0:
                    valid_anchors.append(sentinel_anchor)

        self.dump_to_jsonl(valid_anchors, out_file_path)

    def dump_dictionary_of_counts(self, in_dictionary_path: str)-> None:
        """
        It rewrites the anchor dictionary that contains positions by adding the number of reads alinged to each anchor.

        Parameters
        ----------
        in_dictionary_path : string
            The path to the json object containing the dictionary.
        """
        if not(os.path.exists(in_dictionary_path)):
            print(f"WARNING. Could not find {in_dictionary_path}", file=stderr)
            return

        with open(in_dictionary_path, "r") as f:
            anchors_pos_dict = json.load(f)

        for sentinel in anchors_pos_dict:
            if self.anchor_reads_dict.get(int(sentinel)) == None:
                print(f"something not working here", file=stderr)
                continue
            reads_anchors = self.anchor_reads_dict.get(int(sentinel))

            for i in range(len(anchors_pos_dict[sentinel])):

                # assign to the count the number of reads aligned to the anchor
                anchors_pos_dict[sentinel][i][-1] = len(reads_anchors[i])
        
        #out_dictionary_path = in_dictionary_path[:-5] + ".updated.json"
        self.dump_to_jsonl(anchors_pos_dict, in_dictionary_path)

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

        

    ### PRINT FUNCTION FOR DEBUG ###

    # def get_anchor_string(self, sentinel, index) -> str:
    #     anchor_l = self.sentinel_to_anchor.get(sentinel)
    #     anchor, _ = anchor_l[index]
    #     anchor_str = ""
    #     for node_h in anchor:
    #         orientaiton = "<" if self.graph.get_is_reverse(node_h) else ">"
    #         anchor_str += orientaiton + str(self.graph.get_id(node_h))
    #     return anchor_str
