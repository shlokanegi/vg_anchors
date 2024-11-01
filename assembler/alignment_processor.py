import pickle
import numpy as np
import math

from bdsg.bdsg import PackedGraph
from assembler.constants import READ_P, STRAND_P, START_P, END_P, NODE_P, ORIENT_P, CS_P


class AlignAnchor:

    def __init__(self, packed_graph_path: str, dictionary: dict) -> None:
        #useful initialization objects
        print(packed_graph_path)
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
            # print(f"position: {position} - node {node_id}")
            anchors = self.sentinel_to_anchor.get(node_id)
            if anchors:
                print(f"Found {len(anchors)} anchors: {anchors!r}")
            if anchors:
                for anchor, counter in anchors:
                    # an anchor is a list tuple of a list of node handles 
                    # and a counter set to 0 at the beginning 

                    #locate the position of the sentinel

                    sentinel_p = next(p for p, a  in enumerate(anchor) if self.graph.get_id(a) == node_id)
                    # print(f"Sentinel {node_id} at position {sentinel_p} out of {len(anchor)}")
                    sentinel_orientation = False if self.graph.get_is_reverse(anchor[sentinel_p]) else True
                    concordance_orientation = True if sentinel_orientation == alignment_l[ORIENT_P][position] else False

                    # scan backward to check that the alignment corresponds to the anchor.
                    alignment_matches_anchor = self.verify_path_concordance(position, sentinel_p, concordance_orientation, alignment_l[NODE_P], alignment_l[ORIENT_P], anchor)
                    if alignment_matches_anchor:
                        print(f'Anchor matching {node_id}')
                    
                    # If the node matching is correct, I now need to verify the path is correct

                    # If paths is correct:
                    # 

            walked_length += self.get_length(node_id)
            print(walked_length, flush=True)

            # anchor_size = math.ceil(self.graph.get_length(anchor[0])/2) \
            # + math.floor(self.graph.get_length(anchor[len(anchor) -1 ])/2)

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