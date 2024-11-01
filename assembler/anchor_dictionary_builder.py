from bdsg.bdsg import SnarlDistanceIndex
from bdsg.bdsg import PackedGraph

from assembler.constants import MAX_PATHS_IN_SNARLS, MIN_ANCHOR_LENGTH

import pickle
import math

class SnarlAnchor:

    
    def __init__(self) -> None:
        #useful initialization objects
        self.graph = PackedGraph()
        self.index = SnarlDistanceIndex()
        self.max_num_paths: int = MAX_PATHS_IN_SNARLS
        self.min_anchor_size: int = MIN_ANCHOR_LENGTH

        # important generated_data
        self.leaf_snarls: list = []
        self.sentinel_to_anchor: dict = {}

        # temporary variables to store data between functions
        self.contains_child_snarls: bool = False
        self.keep_path_scan: bool = False
        self.num_nodes: int = 0
        self.temp_nodes: list = []
        self.sentinels: list = []
        self.traversal: list = []
    

    def build_graph(self, packed_graph_path: str,index_path: str )-> None:
        self.graph.deserialize(packed_graph_path)
        self.index.deserialize(index_path)
    

    def check_for_snarl(self,child_net_handle)-> bool:

        if self.index.is_snarl(child_net_handle):
            self.contains_child_snarls = True
        elif self.index.is_node(child_net_handle):
            self.num_nodes += 1
            self.temp_nodes.append(child_net_handle)
        
        # print(f"Children: {self.index.net_handle_as_string(child_net_handle)} has #nodes: {self.num_nodes} and has_snarls is: {self.contains_child_snarls!r}")
        return True
    

    def check_leaf_snarl_iteratee(self, net_handle) -> bool:

        self.contains_child_snarls = False
        self.num_nodes = 0

        # print(f"Visiting {self.index.net_handle_as_string(net_handle)}")

        snarl_children: list = []
        self.index.for_each_child(net_handle, lambda y: snarl_children.append(y) or True)

        self.temp_nodes = []
        for s_c in snarl_children:
            self.index.for_each_child(s_c, self.check_for_snarl)

        if not self.contains_child_snarls: # and self.num_nodes < self.max_num_nodes
            self.leaf_snarls.append(net_handle) # (net_handle, self.temp_nodes)

        return True
    

    # Traverse the whole snalr tree and store the leaf snarls 
    def process_snarls(self) -> None:
        self.index.traverse_decomposition( 
            self.check_leaf_snarl_iteratee,    # snarl_iteratee
            lambda x: True, #  chain_iteratee
            lambda y: True  # node_iteratee
        )
        self.temp_nodes = []
        return
    

    # Returns the list of leaf snalrs (their net_handle)
    def get_leaf_snarls(self)-> list:
        if len(self.leaf_snarls) == 0:
            self.process_snarls()
        return self.leaf_snarls
    

    def get_paths_traversing_snarl(self,snarl_net_handle) -> None:
        start_node_handle, end_node_handle = self.get_snarl_boundaries_handle(snarl_net_handle)

        # print(f"Obtaining steps on start node for snarl {snarl_net_handle}... ", end="", flush=True)
        steps_on_start_node: list = []
        self.graph.for_each_step_on_handle(start_node_handle,lambda y: steps_on_start_node.append(y) or True)
        # print("Done.", flush=True)

        if len(steps_on_start_node) > self.max_num_paths:
            # print(f"Too many paths: found {len(steps_on_start_node)} > {self.max_num_paths} (max paths)", flush=True)
            # for step in steps_on_start_node:
                # path_handle = self.graph.get_path_handle_of_step(step)
                # path_name = self.graph.get_path_name(path_handle)
                # print(path_name, end=" - ")
            # print()
            return []

        # print(f"Obtaining path handles from steps... ", end="", flush=True)
        path_handles: list = []
        for step in steps_on_start_node:
            path_handle = self.graph.get_path_handle_of_step(step)
            path_handles.append(path_handle)
        # print("Done.", flush=True)

        start_is_reverse = self.graph.get_is_reverse(start_node_handle)
        self.sentinels = [self.graph.get_id(end_node_handle), self.graph.get_id(start_node_handle)]\
        if start_is_reverse else [self.graph.get_id(start_node_handle), self.graph.get_id(end_node_handle)]

        snarl_traversals: list = []

        # print(f"Obtaining path in snarls... ", end="", flush=True)
        for path_handle in path_handles:
            self.graph.for_each_step_in_path(path_handle, self.traverse_step_iteratee)
            if self.get_anchor_size(self.traversal) >= self.min_anchor_size:
                # print("anchor inserted", flush=True)
                snarl_traversals.append(self.traversal)
            # else:
                # print(f"anchor of length {self.get_anchor_size(self.traversal)} is too short", flush=True)
            self.traversal = []
        # print("Done.", flush=True)
        return snarl_traversals
    
    def fill_anchor_sentinel_table_single_snarl(self,snarl_net_handle) -> None:
        anchors_list: list = self.get_paths_traversing_snarl(snarl_net_handle)

        for anchor in anchors_list:
            sentinel: int = self.get_sentinel_id(anchor)

            if sentinel not in self.sentinel_to_anchor:
                self.sentinel_to_anchor[sentinel] = [(anchor,[])]
            else:
                insert = True
                for inserted_anchor in self.sentinel_to_anchor[sentinel]:
                    if self.is_equal_anchor(anchor, inserted_anchor):
                        insert = False
                        break
                if insert:
                    self.sentinel_to_anchor[sentinel].append((anchor,[]))
                    # appending a touple containing the anchor and a vector to store the reads it is aligned to

    def fill_anchor_sentinel_table(self) -> None:
        if len(self.leaf_snarls) == 0:
            self.process_snarls()

        for snarl_net_h in self.leaf_snarls:
            self.fill_anchor_sentinel_table_single_snarl(snarl_net_h)


    def get_sentinel_id(self, traversal: list):
        sentinel_handle = traversal[len(traversal)//2]
        return self.graph.get_id(sentinel_handle)
    

    # Returns a tuple containing the start and end boundaries of a snarl
    def get_snarl_boundaries_handle(self, snarl_net_handle) -> tuple:
        start_bound = self.index.get_start_bound(snarl_net_handle)
        end_bound = self.index.get_end_bound(snarl_net_handle)
        return (self.index.get_handle(start_bound,self.graph),self.index.get_handle(end_bound, self.graph))
    

    def traverse_step_iteratee(self, step_handle) -> bool:
        node_handle = self.graph.get_handle_of_step(step_handle)
        node_id = self.graph.get_id(node_handle)

        if not self.keep_path_scan and node_id not in self.sentinels:
            # skipping this node as not in the snalr
            return True
        
        self.traversal.append(node_handle)

        if self.keep_path_scan and node_id in self.sentinels:
            # restore the variable for path scanning to false and return False to stop the iteration
            self.keep_path_scan = False
            return False

        #if keep is false and the node_id is in the sentinels I set it to True.
        #if keep is true and the node_id is not in the sentinels I keep it to True
        self.keep_path_scan = True

        #returning True to keep the iteration going
        return True


    def is_equal_anchor(self, path_1: list, path_2: list) -> bool:
        # if different in length, false
        if len(path_1) != len(path_2):
            return False
        
        #assigning starting postion depending on orientation
        start_orientation_1 = self.graph.get_is_reverse(path_1[0])
        start_orientation_2 = self.graph.get_is_reverse(path_2[0])
        pos_path_1: int = len(path_1) - 1 if start_orientation_1 else 0
        pos_path_2: int = len(path_2) - 1 if start_orientation_2 else 0

        orientation_concordance = start_orientation_1 and start_orientation_2

        for i in range(len(path_1)):
            node_id_1 = self.graph.get_id(path_1[pos_path_1])
            node_id_2 = self.graph.get_id(path_2[pos_path_2])
            if node_id_1 == node_id_2:
                orientation_1 = self.graph.get_is_reverse(path_1[pos_path_1])
                orientation_2 = self.graph.get_is_reverse(path_2[pos_path_2])
                if (orientation_1 and orientation_2) == orientation_concordance:
                    continue
                else: 
                    return False
            else:
                return False
        return True


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
    
    def print_anchors_from_dict(self) -> None:
        for sentinel, anchor_list in self.sentinel_to_anchor.items():
            for anchor, path_l in anchor_list:
                anchor_str = ""
                bandage_nodes_str = ""
                for node_h in anchor:
                    orientaiton = "<" if self.graph.get_is_reverse(node_h) else ">"
                    anchor_str += orientaiton + str(self.graph.get_id(node_h))
                    bandage_nodes_str += "," + str(self.graph.get_id(node_h))
                print(f"Sentinel: {sentinel} ; Anchor : {anchor_str} ; Bandage : {bandage_nodes_str[1:]}")


    def print_traversal(self, traversal: list) -> None:
        for node_handle in traversal:
            #node_handle = self.index.get_handle(node_net_handle, self.graph)
            direction = "<" if self.graph.get_is_reverse(node_handle) == True else ">"
            print(f"{direction}{self.graph.get_id(node_handle)}",end="")
        print()

    def simple_snarl_iteratee(self,snarl_handle) -> bool:
        print("Snarl:", self.index.net_handle_as_string(snarl_handle))
        return True


    def simple_chain_iteratee(self,chain_handle) -> bool:
        print("Chain:", self.index.net_handle_as_string(chain_handle))
        return True


    def simple_node_iteratee(self,node_handle) -> bool:
        print("Node:", self.index.net_handle_as_string(node_handle))
        return True
    
    def print_tree_structure(self) -> None:
        self.index.traverse_decomposition( 
            self.simple_snarl_iteratee,    # snarl_iteratee
            self.simple_chain_iteratee, #  chain_iteratee
            self.simple_node_iteratee  # node_iteratee
        )


    def store_sentinel_dict(self, path : str = "sentinel_to_anchros.pickle") -> None:
        with open(path, 'wb') as f:
            pickle.dump(self.sentinel_to_anchor, f)


    def load_sentinel_dict(self, path : str = "sentinel_to_anchros.pickle") -> None:
        with open(path, 'rb') as f:
            self.sentinel_to_anchor = pickle.load(f)
    
    def get_dict(self) -> dict:
        return self.sentinel_to_anchor