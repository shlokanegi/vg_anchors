import pickle
import numpy as np
import math

from bdsg.bdsg import PackedGraph
from assembler.constants import READ_P, STRAND_P, START_P, END_P, NODE_P, ORIENT_P, CS_P


class AlignAnchor:

    def __init__(self, packed_graph_path: str, dictionary: dict) -> None:
        #useful initialization objects
        self.graph = PackedGraph()
        self.graph.deserialize(packed_graph_path)
        self.sentinel_to_anchor: dict = dictionary

        # temporary variables to store data between functions
    
    def get_length(self, node_id: int) -> int:
        # print(f"getting length of {node_id}")
        node_handle = self.graph.get_handle(node_id)
        # print("done")
        length = self.graph.get_length(node_handle)
        return length
    
    def process_alignment(self, alignment_l):
        # walk through the nodes keeping the position in the node list
        # and the total length of the nodes
        walked_length = 0
        # print(alignment_l[NODE_P])
        # print(alignment_l[ORIENT_P])
        for position, node_id in enumerate(alignment_l[NODE_P]):
            print(f"position: {position} - node {node_id}")
            anchors = self.sentinel_to_anchor.get(node_id)
            # if anchors:
            #     print(f"Found {len(anchors)} anchors: {anchors!r}")
            if anchors:
                for index, anchor_capsule in enumerate(anchors):
                    # an anchor is a list tuple of a list of node handles 
                    # and a counter set to 0 at the beginning 
                    anchor = anchor_capsule[0]
                    #locate the position of the sentinel

                    sentinel_p = next(p for p, a  in enumerate(anchor) if self.graph.get_id(a) == node_id)
                    print(f"Sentinel {node_id} at position {sentinel_p} out of {len(anchor)}")
                    sentinel_orientation = False if self.graph.get_is_reverse(anchor[sentinel_p]) else True
                    concordance_orientation = True if sentinel_orientation == alignment_l[ORIENT_P][position] else False

                    # scan backward to check that the alignment corresponds to the anchor.
                    alignment_matches_anchor = self.verify_path_concordance(position, sentinel_p, concordance_orientation, alignment_l[NODE_P], alignment_l[ORIENT_P], anchor)
                    if alignment_matches_anchor:
                        print(f'Anchor matching {node_id}')
                    
                    # If the node matching is correct, I now need to verify the path is correct
                    anchor_start, anchor_end = self.get_anchor_boundaries(walked_length, anchor, sentinel_p)
                    print(f"Anchor start: {anchor_start} ; end: {anchor_end} (absolute coordinates),{anchor_end-anchor_start} anchor_size {self.get_anchor_size(anchor)}")
                    print(f"al_cs: {alignment_l[CS_P]}, strand: {alignment_l[STRAND_P]}, al_s: {alignment_l[START_P]}, al_end: {alignment_l[END_P]}")

                    x = (anchor_start, anchor_end, alignment_l[CS_P], alignment_l[STRAND_P], alignment_l[START_P], alignment_l[END_P])
                    print(f"I see {len(x)} elements.")
                    is_aligning, read_start, read_end = self.verify_sequence_agreement(*x)
                    # If paths is correct:
                    # I need to append the read info to the anchor hihi.
                    # I need read start and read end of the anchor and the orientation of the read
                    if is_aligning:
                        anchors[index][1].append([alignment_l[READ_P], alignment_l[START_P], read_start, read_end])
                        print(f'Found alignment in read {alignment_l[READ_P]} from {read_start} to {read_end} (length {read_end - read_start}, anchor_size {self.get_anchor_size(anchor)})')
                        # found, no need to check in other anchors
                        break
            walked_length += self.get_length(node_id)
            # print(walked_length, flush=True)

    def get_anchor_boundaries(self, walked_length,  anchor, sentinel_p):
        """
        Input: The total base pairs walked till the sentinel_node, the anchor node_handles list, the position of the sentinel in the anchor list.
        Output: a tuple with the absolute position of the anchor in the path [start,end] start and end are included
        """
        # I need to remember what I chose when building anchor info:
        # if a node has odd size, it is goind to be divided as [l,l+1].
        # so for the boundary nodes, in case they are odds, I will take 
        # the last l + 1 nucleotides of the starting node
        # the first l nucleotides of the last node
        walk_backward = math.ceil(self.graph.get_length(anchor[0])/2) # anchor_length_in_first_boundary_node
        walk_forward = math.floor(self.graph.get_length(anchor[-1])/2) # anchor_length_in_last_boundary_node

        for node_h in anchor[1:sentinel_p]:
            walk_backward += self.graph.get_length(node_h)

        for node_h in anchor[sentinel_p:-1]:
            walk_forward += self.graph.get_length(node_h)
        
        return (walked_length - walk_backward, walked_length + walk_forward)
        

    def verify_path_concordance(self, position, sentinel_p, concordance_orientation, node_array, orientation_array, anchor) -> bool:

        # if they concorde
        if concordance_orientation:
            an_pos = 0
            al_pos = position - sentinel_p 

            if al_pos < 0:
                return False

            while (an_pos < len(anchor) and al_pos < len(node_array)):
                #check node_id is the same
                print(f"al_n: {node_array[al_pos]} | ac_n: {self.graph.get_id(anchor[an_pos])}")
                if node_array[al_pos] != self.graph.get_id(anchor[an_pos]):
                    return False
        
                #check concordance is mantained
                print(f"al_o: {orientation_array[al_pos]} | ac_o: {not self.graph.get_is_reverse(anchor[an_pos])} | conc: {concordance_orientation}")
                if concordance_orientation != (orientation_array[al_pos] == (not self.graph.get_is_reverse(anchor[an_pos]))):
                    return False

                an_pos += 1
                al_pos += 1
        
        else:
            an_pos = 0
            al_pos = position + sentinel_p

            while (an_pos < len(anchor) and al_pos >=0):
                #check node_id is the same
                print(f"al_n: {node_array[al_pos]} | ac_n: {self.graph.get_id(anchor[an_pos])}")
                if node_array[al_pos] != self.graph.get_id(anchor[an_pos]):
                    return False
        
                #check concordance is mantained
                print(f"al_o: {orientation_array[al_pos]} | ac_o: {(not self.graph.get_is_reverse(anchor[an_pos]))} | conc: {concordance_orientation}")
                if concordance_orientation != (orientation_array[al_pos] == (not self.graph.get_is_reverse(anchor[an_pos]))):
                    return False

                an_pos += 1
                al_pos -= 1

        return True
    
# x = (anchor_start, anchor_end, alignment_l[CS_P], alignment_l[STRAND_P], alignment_l[START_P], alignment_l[END_P])
    def verify_sequence_agreement(self, anchor_bp_start, anchor_bp_end, cs_walk, seq_strand, start_in_path, end_in_path):
        # print(seq_strand)
        # If anchor overflows the alingment, it is not valid
        if anchor_bp_end > end_in_path or anchor_bp_start < start_in_path:
            return (False, 0, 0)
        
        walked_in_the_sequence: int = 0 # I need this to keep track of anchor position in the sequence
        walked_in_the_path: int = 0 # I need this to keep track of my walk in the path
        allow_seq_diff: bool = True # I need this to control no variation beteen anchor and sequence is present. Starting with True, setting to False when walking on anchor coordinates

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
                # If I passed on a equal step, it is ok. I set allo_differences to false and go on. But before I check if I have surpassed the end of the anchor. If yes return true.
                if walked_in_the_path >= anchor_bp_end:
                    diff_start = walked_in_the_path - anchor_bp_start
                    diff_end = walked_in_the_path - anchor_bp_end
                    # I add a + 1 in the read_end position because of Shasta requirement that the interval is open at the end. The end id in the sequence is of the first nucleotide after the anchor
                    return (True, walked_in_the_sequence - diff_start, walked_in_the_sequence - diff_end + 1)
                else:
                    allow_seq_diff = False # go to the next step

            # Walking in the anchor section and found a diff
            elif not(allow_seq_diff) and step[0] != ":":
                return (False, 0, 0) 
            
            # I passed the end of the scan and there was no difference
            elif walked_in_the_path >= anchor_bp_end:
                diff_start = walked_in_the_path - anchor_bp_start
                diff_end = walked_in_the_path - anchor_bp_end
                # I add a + 1 in the read_end position because of Shasta requirement that the interval is open at the end. The end id in the sequence is of the first nucleotide after the anchor
                return (True, walked_in_the_sequence - diff_start, walked_in_the_sequence - diff_end)

        return None
    
    def get_anchor_size(self, anchor) -> int:
        # if a node has odd size, it is goind to be divided as [l,l+1].
        # so for the boundary nodes, in case they are odds, I will take 
        # the last l + 1 nucleotides of the starting node
        # the first l nucleotides of the last node
        anchor_size = math.ceil(self.graph.get_length(anchor[0])/2) \
            + math.floor(self.graph.get_length(anchor[len(anchor) -1 ])/2)
        
        for node_handle in anchor[1:len(anchor)-1]:
            anchor_size += self.graph.get_length(node_handle)

        return anchor_size