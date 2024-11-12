# bdsg import
from bdsg.bdsg import SnarlDistanceIndex
from bdsg.bdsg import PackedGraph

# package import
from assembler.constants import (
    MAX_PATHS_IN_SNARLS,
    MIN_ANCHOR_LENGTH,
    MIN_NODES_IN_ANCHOR,
)
from assembler.rev_c import rev_c
from assembler.node import Node

# other imports
import time
import json
from sys import stderr
import pickle


class AnchorDictionary:
    """
    This class produces a Dictionary containing anchors in the pangenome graph.
    Anchors are small paths derived from bubbles in the graph (snarls) and are used by Shasta to
    phase and assembly the reads in a sample.
    A path is a succession of nodes and orientations in a bidirected graph, as described in
    https://github.com/lh3/gfatools/blob/master/doc/rGFA.md.
    The anchor dictionary has as a key a "sentinel" node in the anchor path and as value a list of tuples. In each tuple in the first position there is an anchor having the that node as sentinel and in the second position an empty list that will contain the reads associated to the anchor.
    The anchor is a list of node_handles as in the packedgraph implementation of the bdsg library (https://bdsg.readthedocs.io/en/master/index.html).
    NOTE: THIS IS GOING TO CHANGE AS I WON"T STORE ANYMORE THE HANDLES but only their useful info (node_id, orientation, node_length) and the list of reads associated to the anchors will be stored in another dictionary.

    Path example [handle_####,handle_####,handle####] -> >1>2>3 with [>,< ==  node orientation][node_id]
    Anchors example:>1>2>3      Node with ID == 2 is the sentinel of this path
                    >1>2>3>4    Node with ID == 2 is the sentinel of this path
                    <4<3<2<1    Node with ID == 2 is the sentinel of this path
    """

    def __init__(self) -> None:
        # useful initialization objects
        self.graph = PackedGraph()
        self.index = SnarlDistanceIndex()

        # important generated_data
        self.leaf_snarls: list = []
        self.sentinel_to_anchor: dict = {}
        self.main_path = []

        # temporary variables to store data between functions
        self.contains_child_snarls: bool = False
        self.keep_path_scan: bool = False
        # self.num_nodes: int = 0
        # self.temp_nodes: list = []
        self.snarl_boundaries: tuple = ()
        self.traversal: list = []
        self.verbose = False

    def build(self, packed_graph_path: str, index_path: str) -> None:
        """
        Deserializes the packedGraph and SnarlIndexes generates using vg. Does not return anything

        Parameters
        ----------
        packed_graph_path : string
            Path to packedGraph object (.vg)
        index_path : string
            Path to SnarlIndex object (.dist)
        """
        self.graph.deserialize(packed_graph_path)
        self.index.deserialize(index_path)

    def check_snarl_in_children_iteratee(self, child_net_handle) -> bool:
        """
        It iterates on the children of each snarl child (check ) to verify that the snalr does not contain any other snarl and is therefore a snarl leave. If it seesa a snarl it sets the variable contains_child_snarls as true.
        The return True/ False parameter is used to continue or stop the iteration calling this function. It stops if it finds a snarl else continue to check the handles of the snarl childern.

        Parameters
        ----------
        child_net_handle : object
            net_handle object from a SnarlIndex

        Returns
        -------
        bool
        True if iteration has to continue else False
        """

        if self.index.is_snarl(child_net_handle):
            self.contains_child_snarls = True
            return False
        # elif self.index.is_node(child_net_handle):
        # self.num_nodes += 1
        # self.temp_nodes.append(child_net_handle)

        return True

    def check_leaf_snarl_iteratee(self, net_handle) -> bool:
        """
        This function is called on the snarl tree traversal (process_snarls function) when the pointer is on a snarl net_handle. It verifies if the snarl is a leaf snalr (does not contain inside it another snarl like a matrioska).
        If the snarl is found to be a leaf snarl, it appends its net_handle to a list of valid snarl handles to then process them to generate anchors.
        It returns True to keep the iteration going and do not stop it.

        Parameters
        ----------
        net_handle : object
            net_handle object from a SnarlIndex

        Returns
        -------
        bool
        True as iteration has to continue
        """

        self.contains_child_snarls = False
        # self.num_nodes = 0

        snarl_children: list = []
        self.index.for_each_child(
            net_handle, lambda y: snarl_children.append(y) or True
        )

        # self.temp_nodes = []
        for s_c in snarl_children:
            self.index.for_each_child(s_c, self.check_snarl_in_children_iteratee)

        if not self.contains_child_snarls:
            self.leaf_snarls.append(net_handle)

        return True

    def process_snarls(self) -> None:
        """
        This function traverses the whole Snarl Tree index and stores the leaf snarls into a list for future processing into anchors.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.index.traverse_decomposition(
            self.check_leaf_snarl_iteratee,  # snarl_iteratee
            lambda x: True,  #  chain_iteratee
            lambda y: True,  # node_iteratee
        )
        # self.temp_nodes = []
        return None

    def traverse_step_iteratee(self, step_handle) -> bool:
        """
        This function takes a step_handle in the graph and appends the nodes (of the path associated to that step) that are in the snarl. An anchor is a path in a snalr. The list is stored in the object to be used by another function. The snarl_boundaries are the ids of the nodes that just precede or succede the snarl in the graph structure.

        Parameters
        ----------
        step_handle: obj
            step_handl object of a boundary of a snarl

        Returns
        -------
        True if the iteration on the path has to keep going
            False if the iteration on the path has passed the snarl

        """
        node_handle = self.graph.get_handle_of_step(step_handle)
        node_id = self.graph.get_id(node_handle)

        if not self.keep_path_scan and node_id not in self.snarl_boundaries:
            # skipping this node as not in the snalr
            return True
        self.traversal.append(Node(self.graph.get_id(node_handle), self.graph.get_length(node_handle) ,not(self.graph.get_is_reverse(node_handle))))

        if self.keep_path_scan and node_id in self.snarl_boundaries:
            # restore the variable for path scanning to False for the next path traversal
            # and return False to stop the iteration
            self.keep_path_scan = False
            return False

        # if keep is False and the node_id is in the sentinels I set it to True.
        # if keep is True and the node_id is not in the sentinels I keep it to True
        self.keep_path_scan = True

        # returning True to keep the iteration going
        return True

    def get_paths_traversing_snarl(self, snarl_net_handle) -> list:
        """
        This function takes a snarl net_handle and returns a list of accepted paths traversing it that can be used as anchors.
        A path is accepted if
            - it has strictly > 3 nodes in it. This is arbirtrary to avoid some edge cases in other functions when the anchor has just 2 nodes. These anchors would be flagging a delition / not having an insertion.
            -  the length of the sequence spelled by the node labels is > MIN_ANCHOR_LENGTH == MIN_ANCHOR_LENGTH. See config.py to see the value.

        Parameters
        ----------
        snarl_net_handle: obj
        net_handle object of the snarl

        Returns
        -------
        snarl_traversals: list
        List of accepted paths that can be used as anchors.
        """

        start_node_handle, end_node_handle = self.get_snarl_boundaries_handle(
            snarl_net_handle
        )

        path_handles: list = []  # steps_on_start_node
        self.graph.for_each_step_on_handle(
            start_node_handle,
            lambda y: path_handles.append(self.graph.get_path_handle_of_step(y))
            or True)

        if len(path_handles) > MAX_PATHS_IN_SNARLS:
            return []

        start_is_reverse = self.graph.get_is_reverse(start_node_handle)
        self.snarl_boundaries = (
            (self.graph.get_id(end_node_handle), self.graph.get_id(start_node_handle))
            if start_is_reverse
            else (
                self.graph.get_id(start_node_handle),
                self.graph.get_id(end_node_handle),
            )
        )

        snarl_traversals: list = []

        for path_handle in path_handles:
            self.graph.for_each_step_in_path(path_handle, self.traverse_step_iteratee)
            if (
                self.get_anchor_size(self.traversal) >= MIN_ANCHOR_LENGTH
                and len(self.traversal) >= MIN_NODES_IN_ANCHOR
            ):
                snarl_traversals.append(self.traversal)

            self.traversal = []

        return snarl_traversals

    def fill_anchor_dictionary_single_snarl(self, snarl_net_handle) -> None:
        """
        This function takes a snarl net_handle and fills the sentinel_to_anchor dictionary with the anchors associated to the snarl.

        Parameters
        ----------
        snarl_net_handle: obj
        net_handle object of the snarl

        Returns
        -------
        None
        """
        anchors_list: list = self.get_paths_traversing_snarl(snarl_net_handle)

        for anchor in anchors_list:
            sentinel: int = self.get_sentinel_id(anchor)
            if sentinel == 47:
                print(f"sentinel {sentinel}, anchor: {anchor}")

            if sentinel not in self.sentinel_to_anchor:
                self.sentinel_to_anchor[sentinel] = [anchor]

            else:
                insert = True
                for inserted_anchor in self.sentinel_to_anchor[sentinel]:
                    # verify that the anchor is not already existing in the dictionary
                    if self.is_equal_path(anchor, inserted_anchor):
                        insert = False
                        break
                if insert:
                    self.sentinel_to_anchor[sentinel].append(anchor)

    def fill_anchor_dictionary(self) -> None:
        """
        This function fills the sentinel_to_anchor dictionary with the anchors associated to all the leaf snarls in the graph.

        Returns
        -------
        None
        """
        # in case the leaf snarls were not already computed
        if len(self.leaf_snarls) == 0:
            self.process_snarls()

        for idx, snarl_net_h in enumerate(self.leaf_snarls):
            t0 = time.time()
            self.fill_anchor_dictionary_single_snarl(snarl_net_h)
            print(
                f"Processed snarl {idx}/{len(self.leaf_snarls)} in {time.time()-t0:.2f}",
                file=stderr,
            )

    ### HELPER FUNCTIONS ###

    def get_sentinel_id(self, traversal: list):
        """
        This function takes a path traversing a snarl (candidate anchor) and returns the sentinel associated to it. The sentinel is not dependent on the anchor orientation.

        Parameters
        ----------
        traversal: list
        the list of node handles of the path

        Returns
        -------
        sentinel: int
        the node_id of the sentinel node
        """
        tr_len = len(traversal)
    
        # Determine the index based on the orientation of the first element
        if traversal[0].orientation:
            # If reversed, count from the end
            index = -((tr_len + 1) // 2)
        else:
            # If not reversed, count from the beginning
            index = (tr_len - 1) // 2
    
        return traversal[index].id
    
    # Returns a tuple containing the start and end boundaries of a snarl
    def get_snarl_boundaries_handle(self, snarl_net_handle) -> tuple:
        """
        This function takes a snarl net_handle and returns the boundary nodes of the snarl, i.e. preceding and succeding the snarl. This is used in the candidate anchor generation when traversing the paths to record only the portion of path in the snarl.

        Parameters
        ----------
        snarl_net_handle: obj
        net_handle object of the snarl

        Returns
        -------
        boundary : tuple
        the node_handle of the nodes preceding and succeding the snarl
        """

        start_bound = self.index.get_start_bound(snarl_net_handle)
        end_bound = self.index.get_end_bound(snarl_net_handle)
        boundary = (
            self.index.get_handle(start_bound, self.graph),
            self.index.get_handle(end_bound, self.graph),
        )
        return boundary

    def is_equal_path(self, path_1: list, path_2: list) -> bool:
        """
        This function verifies if two paths in the same snarls are equal, independently of the orientation.

        Parameters
        ----------
        path_1: list
        list of node handle objects of the first path
        path_2: list
        list of node handle objects of the first path

        Returns
        -------
        bool
        True if the paths are equal else False
        """

        # if different in length, don't need to bother
        if len(path_1) != len(path_2):
            return False

        pos_1 = pos_2 = 0
        orientation_concordance = path_1[0].orientation == path_2[0].orientation
        if not orientation_concordance:
            pos_1 = len(path_1) - 1
        
        for _ in range(len(path_1)):
            if path_1[pos_1].id != path_2[pos_2].id or orientation_concordance != (path_1[pos_1].orientation == path_2[pos_2].orientation):
                return False
            pos_1 += 1 if orientation_concordance else -1
            pos_2 += 1
        return True


    def get_anchor_size(self, anchor: list) -> int:
        """
        This function returns the size, in basepairs, of an anchor.

        Parameters
        ----------
        anchor: list
        list of node handle objects of the anchor path

        Returns
        -------
        anchor_size: int
        The length in basepairs of the anchor
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
    
    def steps_path_iteratee(self, step_handle) -> bool:
        """
        This function is applied to the walk in the path and is used to generate a list of steps, defined as list of nodes and their relative position in the path. The position is given by a coordinate that starts as 0 and then is incremented by the length of the node.

        Parameters
        ----------
        step_handle: obj
            step handle object from path

        Returns
        -------
        True to continue the iteration
        """
        # Get the handle (node) for this step
        handle = self.graph.get_handle_of_step(step_handle)
        
        # Get the node ID
        node_id = self.graph.get_id(handle)
        node_length = self.graph.get_length(handle)

        tot_len = 0
        if len(self.main_path) > 0:
            tot_len = self.main_path[-1][-1]

        self.main_path.append((node_id,tot_len + node_length))

        return True
    
    def dump_dictionary(self, out_file_path:str)-> None:
        with open(out_file_path, 'wb') as out_f:
            pickle.dump(self.sentinel_to_anchor, out_f)

    
    def generate_positioned_dictionary(self, graph_path_name: str, out_file_path:str) -> None:
        """
        This function generates a dictionary mirroring the sentinel_to_anchor dictionary but instead it populates the anchors with node_ids associated to the anchor and the minimum position of the nodes it is made of. 

        Parameters
        ----------
        graph_path_name: string
            The name of the path used as referencing for the position (I use CHM13)
        out_file_path: string
            The path of the output json where to store the dictionary

        Returns
        -------
        None
        """
        graph_path_name = "CHM13#chr20:149948-250000"

        if self.graph.has_path(graph_path_name):
            path_handle = self.graph.get_path_handle(graph_path_name)
        
        # generate main path
        self.graph.for_each_step_in_path(path_handle, self.steps_path_iteratee)
        size_dict = dict()

        for node_id, pos in self.main_path:
            size_dict[node_id] = pos

        anchor_pos_dict = dict()
        for sentinel, anchor_list in self.sentinel_to_anchor.items():
            for anchor in anchor_list:
                min_pos = max( size_dict.get(x.id,-1) for x in anchor)
                if anchor_pos_dict.get(sentinel) == None:
                    anchor_pos_dict[sentinel] = []
                anchor_pos_dict[sentinel].append(([x.id for x in anchor],min_pos,0))
        
        with open(out_file_path, "w", encoding="utf-8") as f:
            json.dump(anchor_pos_dict, f, ensure_ascii=False)


    ### PRINTING FUNCTIONS FOR DEBUG - VISUALIZATION ###

    def print_anchors_from_dict(self, file_name) -> None:
        with open(file_name, "w") as f:
            for sentinel, anchor_list in self.sentinel_to_anchor.items():
                for anchor in anchor_list:
                    anchor_str = ""
                    bandage_nodes_str = ""
                    for node in anchor:
                        orientaiton = ">" if node.orientation else "<"
                        anchor_str += orientaiton + str(node.id)
                        bandage_nodes_str += "," + str(node.id)
                    print(
                        f"Sentinel: {sentinel} ; Anchor : {anchor_str} ; Bandage : {bandage_nodes_str[1:]}",
                        file=f,
                    )

    def print_traversal(self, traversal: list) -> None:
        for node in traversal:
            # node_handle = self.index.get_handle(node_net_handle, self.graph)
            direction = ">" if node.orientation == True else "<"
            print(f"{direction}{node.id}", end="")
        print()

    def print_tree_snarl_iteratee(self, snarl_handle) -> bool:
        print("Snarl:", self.index.net_handle_as_string(snarl_handle), flush=True)
        return True

    def print_tree_chain_iteratee(self, chain_handle) -> bool:
        print("Chain:", self.index.net_handle_as_string(chain_handle), flush=True)
        return True

    def print_tree_node_iteratee(self, node_handle) -> bool:
        print("Node:", self.index.net_handle_as_string(node_handle), flush=True)
        return True

    def print_tree_structure(self) -> None:
        print("Printing tree structure now", flush=True)
        self.index.traverse_decomposition(
            self.print_tree_snarl_iteratee,  # snarl_iteratee
            self.print_tree_chain_iteratee,  #  chain_iteratee
            self.print_tree_node_iteratee,  # node_iteratee
        )

    def print_sentinels_for_bandage(self, file) -> None:
        sentinel_nodes_set = set()
        for _, anchor_list in self.sentinel_to_anchor.items():
            for anchor in anchor_list:
                for node_h in anchor:
                    sentinel_nodes_set.add(node_h.id)

        with open(file, "w") as out_f:
            print("Node,color", file=out_f)
            for node in sentinel_nodes_set:
                print(f"{node},#FF0000", file=out_f)

    def get_dict(self) -> dict:
        return self.sentinel_to_anchor

    # def get_anchor_seq(self, anchor) -> str:
    #     anchor_seq = ""

    #     orientation = anchor[0].orientation
    #     anchor_c = anchor[:] if orientation else anchor[::-1]

    #     p_start_node = anchor_c[0].length // 2
    #     p_end_node = anchor_c[-1].length // 2
    #     anchor_size = p_start_node + p_end_node

    #     node_f_seq = self.graph.get_sequence(anchor_c[0])
    #     node_f_seq = rev_c(node_f_seq) if reversed_orientation else node_f_seq
    #     anchor_seq += node_f_seq[self.graph.get_length(anchor_c[0]) - p_start_node :]

    #     for node_handle in anchor_c[1:-1]:
    #         anchor_size += self.graph.get_length(node_handle)
    #         s = self.graph.get_sequence(node_handle)
    #         anchor_seq += rev_c(s) if reversed_orientation else s

    #     node_f_seq = self.graph.get_sequence(anchor[-1])
    #     node_f_seq = (
    #         rev_c(node_f_seq) if self.graph.get_is_reverse(anchor_c[-1]) else node_f_seq
    #     )
    #     anchor_seq += node_f_seq[:p_end_node]

    #     print(f"Anchor seq size == {len(anchor_seq)}, anchor_size == {anchor_size}")
    #     return anchor_seq

    def print_dict_sizes(self, out_f) -> None:
        with open(out_f, "w") as f:
            for sentinel, anchor_list in self.sentinel_to_anchor.items():
                for anchor in anchor_list:
                    # anchor_s = self.get_anchor_seq(anchor)
                    print(
                        f"{sentinel},{self.get_anchor_size(anchor)}", 
                        file=f,
                    )
            #{anchor_s},{rev_c(anchor_s)}

    # Returns the list of leaf snalrs (their net_handle)

    # def get_leaf_snarls(self) -> list:
    #     if len(self.leaf_snarls) == 0:
    #         self.process_snarls()
    #     return self.leaf_snarls

    # Serializing the dictionary using pickle - can't be use with node hanles

    # def store_sentinel_dict(self, path: str = "sentinel_to_anchros.pickle") -> None:
    #     with open(path, "wb") as f:
    #         pickle.dump(self.sentinel_to_anchor, f)

    # def load_sentinel_dict(self, path: str = "sentinel_to_anchros.pickle") -> None:
    #     with open(path, "rb") as f:
    #         self.sentinel_to_anchor = pickle.load(f)
