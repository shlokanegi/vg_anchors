# bdsg import
# import json
import time
from bdsg.bdsg import SnarlDistanceIndex
from bdsg.bdsg import PackedGraph

# package import
from assembler.constants import (
    # MAX_PATHS_IN_SNARLS,
    MIN_ANCHOR_LENGTH,
    MIN_NODES_IN_ANCHOR,
    FORWARD_DICTIONARY,
    REVERSE_DICTIONARY,
    PEEK_SIZE,
    END_NODE_POS,
    SNARL_ID_POS,
)
from assembler.rev_c import rev_c
from assembler.node import Node
from assembler.anchor import Anchor

# other imports
import time
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
        self.snarl_boundaries: list = [dict(), dict()]
        
        # temporary variables to store data between functions
        self.contains_child_snarls: bool = False
        self.keep_path_scan: bool = True
        self.path_orientation: bool = True
        self.current_snarl_start: int = 0
        self.peek_orientations = []
        self.count_in_path: bool = True
        self.num_usable_bubbles = 0
        self.next_handle_expand_boundary = None
        self.anchor_length_occupied = 0

        self.current_anchor: Anchor = Anchor()
        self.curr_path_name = ""
        self.verbose = False
        self.ref_path_name = "CHM13".casefold()
        self.path_names = []

        # variables used for debugging
        self.used_bubbles = dict()

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
        t0=time.time()
        self.graph.deserialize(packed_graph_path)
        self.index.deserialize(index_path)
        print(f"Graph files loaded in {time.time()-t0:.2f}", file=stderr)

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

        snarl_children: list = []
        self.index.for_each_child(
            net_handle, lambda y: snarl_children.append(y) or True
        )

        for snarl_child in snarl_children:
            self.index.for_each_child(snarl_child, self.check_snarl_in_children_iteratee)

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
        return None

    def get_path_orientation_iteratee(self, step_handle) -> bool:
        node_handle = self.graph.get_handle_of_step(step_handle)
        orientation = not self.graph.get_is_reverse(node_handle)
        # if self.graph.get_is_reverse(node_handle):
        #     orientation = False
        self.peek_orientations.append(orientation)
        # if len(self.peek_orientations) >= PEEK_SIZE:
        #     num_forward = self.peek_orientations.count(True)
        #     self.path_orientation = (
        #         FORWARD_DICTIONARY if num_forward >= (PEEK_SIZE // 2) else REVERSE_DICTIONARY
        #     )
        #     self.peek_orientations = []
        #     return False
        return True


    def traverse_step_iteratee(self, step_handle) -> bool:
        """
        This function takes a step_handle in the graph and appends the nodes (of the path associated to that step) that are in the snarl. An anchor is a path in a snarl. 
        The list is stored in the object to be used by another function. The snarl_boundaries are the ids of the nodes that just precede or succeed the snarl in the graph structure.
        If an anchor is found, the sentinel_to_anchor dictionary is also updated.

        Parameters
        ----------
        step_handle: obj
            step_handle object of a boundary of a snarl

        Returns
        -------
        True if the iteration on the path has to keep going
            False if the iteration on the path has passed the snarl

        """
        node_handle = self.graph.get_handle_of_step(step_handle)
        node_id = self.graph.get_id(node_handle)

        print(f"In current path, traversing node id {node_id}", end="\n")

        if (
            not self.keep_path_scan
            and self.snarl_boundaries[self.path_orientation][self.current_snarl_start][END_NODE_POS] != node_id
            and node_id in self.snarl_boundaries[self.path_orientation][self.current_snarl_start][2]
        ):
            self.current_anchor.add(
                Node(
                    self.graph.get_id(node_handle),
                    self.graph.get_length(node_handle),
                    not (self.graph.get_is_reverse(node_handle)),
                )
            )
            print(f"Adding node {self.graph.get_id(node_handle)} to anchor gets {self.current_anchor!r}", end="\n")

            return True

        elif (
            not self.keep_path_scan
        ):
            if self.snarl_boundaries[self.path_orientation][self.current_snarl_start][END_NODE_POS] == node_id:
                self.current_anchor.add(
                    Node(
                        self.graph.get_id(node_handle),
                        self.graph.get_length(node_handle),
                        not (self.graph.get_is_reverse(node_handle)),
                    )
                )
                
                self.current_anchor.compute_snarl_boundary()
                self.current_anchor.bp_occupied_start_node = (0 if self.current_anchor.snarl_start_node.length == 1 else 1)
                self.current_anchor.bp_occupied_end_node = (0 if self.current_anchor.snarl_end_node.length == 1 else 1)
                self.current_anchor.compute_bp_length()
                
                self.current_anchor.compute_sentinel_bp_length()
                
                if (
                    len(self.current_anchor) >= MIN_NODES_IN_ANCHOR
                ):
                    sentinel: int = self.current_anchor.get_sentinel_id()
                    # print(f"Sentinel is {sentinel}", file=stderr)
                    if sentinel not in self.sentinel_to_anchor:
                        # print(f"Added new anchor at sentinel {sentinel}", file=stderr)
                        self.current_anchor.add_reference_path(self.curr_path_name)
                        self.sentinel_to_anchor[sentinel] = [self.current_anchor]
                        print(f"Final anchor is {self.current_anchor!r} whose sentinal is {sentinel} and length {self.current_anchor.basepairlength}")

                    else:
                        insert = True
                        for id, inserted_anchor in enumerate(self.sentinel_to_anchor[sentinel]):
                            # verify that the anchor is not already existing in the dictionary
                            if self.current_anchor == inserted_anchor:  # if current anchor is same as the one already in the list at index 'id', then just update the path variable of the anchor
                                self.sentinel_to_anchor[sentinel][id].add_reference_path(self.curr_path_name)
                                print(f"Final anchor is {self.current_anchor!r} whose sentinal is {sentinel} and length {self.current_anchor.basepairlength}")
                                insert = False
                        if insert:
                            # but, if current anchor is not already in the list, then add it to the list and also update path variable
                            self.current_anchor.add_reference_path(self.curr_path_name)
                            self.sentinel_to_anchor[sentinel].append(self.current_anchor)
                            print(f"Final anchor is {self.current_anchor!r} whose sentinal is {sentinel} and length {self.current_anchor.basepairlength}")
            self.current_anchor = Anchor()
            self.keep_path_scan = True

        if (
            self.keep_path_scan
            and self.snarl_boundaries[self.path_orientation].get(node_id) != None
        ):

            self.current_snarl_start = node_id
            self.current_anchor.add(
                Node(
                    self.graph.get_id(node_handle),
                    self.graph.get_length(node_handle),
                    not (self.graph.get_is_reverse(node_handle)),
                )
            )
            print(f"Adding node {self.graph.get_id(node_handle)} to anchor. Corresponding boundary node is {self.snarl_boundaries[self.path_orientation][self.current_snarl_start][END_NODE_POS]}", end="\n")
            self.current_anchor.add_snarl_id(
                self.snarl_boundaries[self.path_orientation][node_id][SNARL_ID_POS]
            )
            self.used_bubbles[self.snarl_boundaries[self.path_orientation][node_id][SNARL_ID_POS]] = True
            self.keep_path_scan = False
            return True

        # returning True to keep the iteration going
        return True


    def get_edge_snarl(self, snarl_net_handle, extend=False) -> None:
        """
        This function takes a snarl_net_handle (from a list of leaf snarls), computes their boundary nodes along with nodes inside it. It populates FORWARD and REVERSE snarl 
        dictionaries.

        Parameters
        ----------
        snarl_net_handle: obj
            net_handle object of the snarl

        Returns
        -------
        None
        """
        self.num_usable_bubbles += 1
        
        if extend:
            start_node_handle, end_node_handle, nodes_inside = self.get_snarl_boundaries_extend(snarl_net_handle)
        else:
            start_node_handle, end_node_handle, nodes_inside = self.get_snarl_boundaries_handle(snarl_net_handle)

        snarl_boundary = (
            (self.graph.get_id(end_node_handle), self.graph.get_id(start_node_handle))
            if self.graph.get_id(end_node_handle) < self.graph.get_id(start_node_handle)
            else (
                self.graph.get_id(start_node_handle),
                self.graph.get_id(end_node_handle),
            )
        )
        
        self.snarl_boundaries[FORWARD_DICTIONARY][snarl_boundary[0]] = (
            snarl_boundary[1],
            self.num_usable_bubbles,    # SNARL ID
            nodes_inside
        )
        self.snarl_boundaries[REVERSE_DICTIONARY][snarl_boundary[1]] = (
            snarl_boundary[0],
            self.num_usable_bubbles,
            nodes_inside
        )
        return

    def print_anchor_boundaries_stderr(self):
        for el in self.snarl_boundaries[FORWARD_DICTIONARY]:
            print(
                f"{el},{self.snarl_boundaries[FORWARD_DICTIONARY][el][END_NODE_POS]}", file=stderr
            )

    def print_anchor_boundaries_dict(self, file_path):
        print(f"Printing to {file_path}.forward_dict.csv")
        with open(f"{file_path}.forward_dict.csv", "w") as f:
            for el in self.snarl_boundaries[FORWARD_DICTIONARY]:
                print(
                    f"{el},{self.snarl_boundaries[FORWARD_DICTIONARY][el][END_NODE_POS]},{self.snarl_boundaries[FORWARD_DICTIONARY][el][2]}", file=f
                )
        print(f"Printing to {file_path}.reverse_dict.csv")
        with open(f"{file_path}.reverse_dict.csv", "w") as f:
            for el in self.snarl_boundaries[REVERSE_DICTIONARY]:
                print(
                    f"{el},{self.snarl_boundaries[REVERSE_DICTIONARY][el][END_NODE_POS]},{self.snarl_boundaries[REVERSE_DICTIONARY][el][2]}", file=f
                )

    def collect_path_handles(self, step_handle):
        path_handle = self.graph.get_path_handle_of_step(step_handle)
        self.path_names.append(self.graph.get_path_name(path_handle))  # self.graph.get_path_name()
        return True

    def get_snalrs_from_paths(self) -> None:
        """
        This function takes a leaf snarl net_handle and fills the sentinel_to_anchor dictionary with the anchors associated to the snarl.

        Parameters
        ----------
        snarl_net_handle: obj
        net_handle object of the snarl

        Returns
        -------
        None
        """

        #collecting path handles to scan in the graph
        for node in self.snarl_boundaries[0]:
            self.graph.for_each_step_on_handle(self.graph.get_handle(node), self.collect_path_handles)

        print(f"TOT PATHS COLLECTED: {len(self.path_names)}")
        self.path_names = sorted(list(set(self.path_names)))
        #scan path handles to obtain the alleles in the snarls.
        print(f"Ready to process {len(self.path_names)} paths...", end = ' ')
        t_0 = time.time()
        for path_name in self.path_names:
            for path_orientation in [REVERSE_DICTIONARY, FORWARD_DICTIONARY]:

                path_handle = self.graph.get_path_handle(path_name)
                self.curr_path_name = path_name
                print(f"Currently processing path {self.curr_path_name}...", end="\n")
                self.current_snarl_start = -1
                self.keep_path_scan = True
                self.count_in_path = True
                self.current_anchor = Anchor()
                self.peek_orientations = []
                self.path_orientation=path_orientation

                print(f"With path_orientation {self.path_orientation}", end="\n")
                self.graph.for_each_step_in_path(path_handle, self.traverse_step_iteratee)

            print(f"done in {time.time()-t_0}")


    def generate_anchors_boundaries(self, extend=False):
        """
        This function sorts leaf snarl handle list based on snarl orientation, so that all snarl handles are in ascending order of occurrence.
        It also updates the snarl dictionaries
        """
        
        ## CHECK if snarls were found in forward direction or reverse
        snarl1_snarl_net_handle, snarl2_snarl_net_handle = self.leaf_snarls[0], self.leaf_snarls[1]

        # snarl1 start node
        snarl1_start_bound_net_handle = self.index.get_start_bound(snarl1_snarl_net_handle)
        snarl1_start_bound_handle = self.index.get_handle(snarl1_start_bound_net_handle, self.graph)
        snarl1_start_node = self.graph.get_id(snarl1_start_bound_handle)
        # snarl2 start node
        snarl2_start_bound_net_handle = self.index.get_start_bound(snarl2_snarl_net_handle)
        snarl2_start_bound_handle = self.index.get_handle(snarl2_start_bound_net_handle, self.graph)
        snarl2_start_node = self.graph.get_id(snarl2_start_bound_handle)

        if snarl1_start_node > snarl2_start_node:
            self.leaf_snarls = self.leaf_snarls[::-1]

        for _, snarl_net_handle in enumerate(self.leaf_snarls):
            self.get_edge_snarl(snarl_net_handle, extend)


    def fill_anchor_dictionary(self, extend = False) -> None:
        """
        This function fills the sentinel_to_anchor dictionary with the anchors associated to all the leaf snarls in the graph.

        Returns
        -------
        None
        """
        # in case the leaf snarls were not already computed
        t0 = time.time()
        if len(self.leaf_snarls) == 0:
            self.process_snarls()
        print(
            f"Leaf Snarls Computed in {time.time()-t0:.2f}",
            file=stderr,
        )

        t1 = time.time()
        self.generate_anchors_boundaries(extend)
        print(
            f"Snarl Boundaries computed in {time.time()-t1:.2f}",
            file=stderr,
        )

        t2 = time.time()
        self.get_snalrs_from_paths()

        print(
            f"Snarl dictionary computed in {time.time()-t2:.2f}. Total time: {time.time()-t0:.2f}.",
            file=stderr,
        )
        num_used_bubbles = len([_ for _ in self.used_bubbles if self.used_bubbles[_] == 1])
        print(f"Num used bubbles: {num_used_bubbles} ; Num usable bubbles: {self.num_usable_bubbles} ; ratio {num_used_bubbles/self.num_usable_bubbles}.",file=stderr)

    ### HELPER FUNCTIONS ###

    def get_nodes_in_snarl(self, snarl_net_handle) -> list:
        nodes_inside = []

        self.index.traverse_decomposition_helper(
            snarl_net_handle,  
            snarl_iteratee=lambda s: True,  # Ignore snarls
            chain_iteratee=lambda c: True,  # Ignore chains
            node_iteratee=lambda n: nodes_inside.append(self.graph.get_id(self.index.get_handle(n, self.graph))) or True
        )

        return nodes_inside  # Return the collected node IDs

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

        start_bound_net_handle = self.index.get_start_bound(snarl_net_handle)
        end_bound_net_handle = self.index.get_end_bound(snarl_net_handle)

        start_bound_handle = self.index.get_handle(start_bound_net_handle, self.graph)
        end_bound_handle = self.index.get_handle(end_bound_net_handle, self.graph)
        nodes_inside_snarl = self.get_nodes_in_snarl(snarl_net_handle)

        return (start_bound_handle,end_bound_handle,nodes_inside_snarl)


    # def get_snarl_boundaries_extend(self, snarl_net_handle) -> tuple:
    #     """
    #     This function takes a snarl net_handle and returns the boundary nodes of the snarl, i.e. preceding and succeding the snarl, allowing for boundary extension in case the boundary nodes length is not enough to generate an anchor. 
    #     This is used in the candidate anchor generation when traversing the paths to record only the portion of path in the snarl.

    #     Parameters
    #     ----------
    #     snarl_net_handle: obj
    #     net_handle object of the snarl

    #     Returns
    #     -------
    #     boundary : tuple
    #     the node_handle of the nodes preceding and succeding the snarl
    #     """

    #     start_bound_net_handle = self.index.get_start_bound(snarl_net_handle)
    #     end_bound_net_handle = self.index.get_end_bound(snarl_net_handle)
    #     nodes_inside_snarl = self.get_nodes_in_snarl(snarl_net_handle)

    #     start_bound_handle = self.index.get_handle(start_bound_net_handle, self.graph)
    #     end_bound_handle = self.index.get_handle(end_bound_net_handle, self.graph)

    #     #TODO: Switch start and end node IDs if end is bigger than start node

    #     print(f"In snarl with boundary ({self.graph.get_id(start_bound_handle)},{self.graph.get_id(end_bound_handle)}).", flush=True)

    #     #1) Verify that the anchor can be built
    #     boundary_nodes_length = self.graph.get_length(start_bound_handle) // 2 + self.graph.get_length(end_bound_handle) // 2
    #     # if the boundary nodes are long enough, return them
    #     if boundary_nodes_length >= MIN_ANCHOR_LENGTH:
    #         return (start_bound_handle,end_bound_handle,nodes_inside_snarl) 
        
    #     # else, try to extend the boundary to nearby nodes, if and only if the node degree is 1.
    #     else:
    #         print(f"Running extention for snarl, since boundary nodes lengths summed up to {boundary_nodes_length}")
            
    #         # Check if node degree of start and/or end node is 1. Try to complete MIN_ANCHOR_LENGTH using primary boundary nodes.
    #         go_left = True
    #         start_degree = self.graph.get_degree(start_bound_handle, go_left)
    #         end_degree = self.graph.get_degree(end_bound_handle, not go_left)
    #         self.anchor_length_occupied = (self.graph.get_length(end_bound_handle) // 2) + (self.graph.get_length(start_bound_handle) // 2)

    #         total = self.anchor_length_occupied \
    #             + ((self.graph.get_length(start_bound_handle) - self.graph.get_length(start_bound_handle) // 2) if start_degree==1 else 0) \
    #             + ((self.graph.get_length(end_bound_handle) - self.graph.get_length(end_bound_handle) // 2) if end_degree==1 else 0)
            
    #         if total >= MIN_ANCHOR_LENGTH:
    #             #1) Find longer node
    #             if self.graph.get_length(start_bound_handle) > self.graph.get_length(end_bound_handle):
                         
            
    #         # First expand boundary on the longest node direction and try to reach MIN_ANCHOR_LENGTH
    #         current_handle = start_bound_handle if go_left else end_bound_handle
    #         computed_handle, nodes_inside_snarl_extended = self.expand_bounary(current_handle, go_left, nodes_inside_snarl)

    #         if self.anchor_length_occupied >= MIN_ANCHOR_LENGTH:
    #             boundary = (computed_handle, end_bound_handle, nodes_inside_snarl) if go_left else (start_bound_handle, computed_handle, nodes_inside_snarl)
    #             return boundary
            
    #         #If not enough, try the other side
    #         go_left = not(go_left)

    #         current_handle = start_bound_handle if go_left else end_bound_handle
    #         other_computed_handle, nodes_inside_snarl_extended = self.expand_bounary(current_handle, go_left, nodes_inside_snarl)
    #         boundary = (other_computed_handle, computed_handle, nodes_inside_snarl) if go_left else (computed_handle, other_computed_handle, nodes_inside_snarl)
    #     return boundary
    
    # def expand_bounary(self, current_handle, go_left_bool, nodes_inside_snarl):
    #     while self.anchor_length_occupied < MIN_ANCHOR_LENGTH:
    #         current_handle_id = self.graph.get_id(current_handle)
    #         print(f" Seeing {current_handle_id}", flush=True)
    #         # if current handle is present in either the forward or reverse dict, i.e. node degree is > 1
    #         if self.snarl_boundaries[FORWARD_DICTIONARY].get(current_handle_id) != None or self.snarl_boundaries[REVERSE_DICTIONARY].get(current_handle_id) != None:
    #             print(f"{current_handle_id} that is in the dictionary. Stopping", flush=True)
    #             break
    #         degree = self.graph.get_degree(current_handle, go_left_bool)

    #         if degree == 1:
    #             self.graph.follow_edges(current_handle, go_left_bool, self.next_handle_iteratee)
    #             if self.next_handle_expand_boundary is None or self.next_handle_expand_boundary == current_handle:
    #                 break
    #             else:
    #                 self.anchor_length_occupied += self.graph.get_length(current_handle) - (self.graph.get_length(current_handle) // 2) + (self.graph.get_length(self.next_handle_expand_boundary) // 2)
    #                 current_handle = self.next_handle_expand_boundary
    #                 nodes_inside_snarl.append(self.graph.get_id(current_handle))

    #         else:
    #             break
    #     return current_handle, nodes_inside_snarl



    def get_snarl_boundaries_extend(self, snarl_net_handle) -> tuple:
        """
        This function takes a snarl net_handle and returns the boundary nodes of the snarl, i.e. preceding and succeding the snarl, allowing for boundary extension in case the boundary nodes length is not enough to generate an anchor. 
        This is used in the candidate anchor generation when traversing the paths to record only the portion of path in the snarl.

        Parameters
        ----------
        snarl_net_handle: obj
        net_handle object of the snarl

        Returns
        -------
        boundary : tuple
        the node_handle of the nodes preceding and succeding the snarl
        """

        start_bound_net_handle = self.index.get_start_bound(snarl_net_handle)
        end_bound_net_handle = self.index.get_end_bound(snarl_net_handle)
        nodes_inside_snarl = self.get_nodes_in_snarl(snarl_net_handle)

        start_bound_handle = self.index.get_handle(start_bound_net_handle, self.graph)
        end_bound_handle = self.index.get_handle(end_bound_net_handle, self.graph)

        print(f"In snarl with boundary ({self.graph.get_id(start_bound_handle)},{self.graph.get_id(end_bound_handle)}).", flush=True)


        #TODO: COMPLETE
        #1) Verify that the anchor can be built
        boundary_nodes_length = self.graph.get_length(start_bound_handle) // 2 + self.graph.get_length(end_bound_handle) // 2
        
        # if the boundary nodes are long enough, return them
        if boundary_nodes_length >= MIN_ANCHOR_LENGTH:
            return (start_bound_handle,end_bound_handle,nodes_inside_snarl) 
        
        # else, try to extend the boundary to nearby nodes, if and only if the node degree is 1.
        else:
            print(f"Running extention for snarl, since boundary nodes lengths summed up to {boundary_nodes_length}")
            go_left = False
            #1) start expansion by the longer node
            if self.graph.get_length(start_bound_handle) > self.graph.get_length(end_bound_handle):
                go_left = True
            
            # self.anchor_length_occupied = self.graph.get_length(end_bound_handle) // 2 if go_left else self.graph.get_length(start_bound_handle) // 2
            self.anchor_length_occupied = (self.graph.get_length(end_bound_handle) // 2) + (self.graph.get_length(start_bound_handle) // 2)
            
            # First expand boundary on the longest node direction and try to reach MIN_ANCHOR_LENGTH
            current_handle = start_bound_handle if go_left else end_bound_handle
            computed_handle, nodes_inside_snarl_extended = self.expand_bounary(current_handle, go_left, nodes_inside_snarl)

            if self.anchor_length_occupied >= MIN_ANCHOR_LENGTH:
                boundary = (computed_handle, end_bound_handle, nodes_inside_snarl_extended) if go_left else (start_bound_handle, computed_handle, nodes_inside_snarl)
                return boundary
            
            #If not enough, try the other side
            go_left = not(go_left)

            current_handle = start_bound_handle if go_left else end_bound_handle
            other_computed_handle, nodes_inside_snarl_extended = self.expand_bounary(current_handle, go_left, nodes_inside_snarl)
            boundary = (other_computed_handle, computed_handle, nodes_inside_snarl_extended) if go_left else (computed_handle, other_computed_handle, nodes_inside_snarl)
        return boundary
    
    def expand_bounary(self, current_handle, go_left_bool, nodes_inside_snarl):
        while self.anchor_length_occupied < MIN_ANCHOR_LENGTH:
            current_handle_id = self.graph.get_id(current_handle)
            print(f" Seeing {current_handle_id}", flush=True)
            # if current handle is present in either the forward or reverse dict
            if self.snarl_boundaries[FORWARD_DICTIONARY].get(current_handle_id) != None or self.snarl_boundaries[REVERSE_DICTIONARY].get(current_handle_id) != None:
                print(f"{current_handle_id} that is in the dictionary. Stopping", flush=True)
                break
            degree = self.graph.get_degree(current_handle, go_left_bool)

            if degree == 1:
                print(f"inside extension, current 1-degree node being checked: {current_handle_id}")
                self.graph.follow_edges(current_handle, go_left_bool, self.next_handle_iteratee)
                if self.next_handle_expand_boundary is None or self.next_handle_expand_boundary == current_handle:
                    break
                else:
                    self.anchor_length_occupied += self.graph.get_length(current_handle) - (self.graph.get_length(current_handle) // 2) + (self.graph.get_length(self.next_handle_expand_boundary) // 2)
                    nodes_inside_snarl.append(self.graph.get_id(current_handle))
                    # nodes_inside_snarl.append(self.graph.get_id(self.next_handle_expand_boundary))
                    current_handle = self.next_handle_expand_boundary

            else:
                break
        return current_handle, nodes_inside_snarl
    

    def next_handle_iteratee(self, next_boundary):
        self.next_handle_expand_boundary = next_boundary
        # returning False as there is just 1 node connected when the degree is 1.
        return False


    def steps_path_iteratee(self, step_handle) -> bool:
        """
        This function is applied to the walk in the path and is used to generate a list of steps, defined as list of nodes and their relative position in the path. 
        The position is given by a coordinate that starts as 0 and then is incremented by the length of the node.

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
            tot_len = self.main_path[-1][-1]       # updates tot_len by the walk till the last node

        self.main_path.append((node_id, tot_len + node_length))

        return True


    def dump_dictionary(self, out_file_path: str) -> None:
        with open(out_file_path, "wb") as out_f:
            pickle.dump(self.sentinel_to_anchor, out_f)


    def add_positions_to_anchors(self, graph_path_name: str = "") -> None:
        """
        This function populates the anchors with their (average) position in the CHM13 path

        Parameters
        ----------
        graph_path_name: string
            The name of the path used as referencing for the position (I use CHM13)

        Returns
        -------
        None
        """

        graph_path_name = ""
        # path_names = self.get_path_names()

        for path in self.path_names:
            if self.ref_path_name in path.casefold():
                graph_path_name = path

        if self.graph.has_path(graph_path_name):
            print(f"Found path {graph_path_name}", file=stderr)
            path_handle = self.graph.get_path_handle(graph_path_name)
        else:
            print(f"WARNING: Could not find CHM13 path in graph", file=stderr)
            return

        # generate main path (node_id, length_in_chm13)
        self.graph.for_each_step_in_path(path_handle, self.steps_path_iteratee)
        size_dict = dict()

        for node_id, pos in self.main_path:
            size_dict[node_id] = pos

        for _, anchor_list in self.sentinel_to_anchor.items():
            for anchor in anchor_list:
                start_node_pos = size_dict.get(anchor[0].id, -1)
                end_node_pos = size_dict.get(anchor[-1].id, -1)
                if start_node_pos > 0 and end_node_pos > 0:
                    anchor.genomic_position = (start_node_pos + end_node_pos) // 2
                else:
                    pos = max(size_dict.get(x.id, -1) for x in anchor)
                    anchor.genomic_position = max(pos, 0)
                anchor.chromosome = graph_path_name

    ### PRINTING FUNCTIONS FOR DEBUG - VISUALIZATION ###

    def print_anchors_from_dict(self, file_name) -> None:
        with open(file_name, "w") as f:
            for sentinel, anchor_list in self.sentinel_to_anchor.items():
                for anchor in anchor_list:
                    print(
                        f"Sentinel: {sentinel} ; Anchor : {anchor!r} ; Bandage : {anchor.bandage_representation()}",
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

    def get_path_names(self):
        path_names = []

        def collect_name(path_handle):
            path_names.append(self.graph.get_path_name(path_handle))
            return True

        self.graph.for_each_path_handle(collect_name)
        return path_names

    def print_dict_sizes(self, out_f) -> None:
        with open(out_f, "w") as f:
            print(f"Sentinel_node\tsnarl_id\tAnchor_length\tAnchor_pos_in_ref_path\tAnchor_path\tAnchor_nodes_copypaste_bandage\tPaths_associated_with_anchor",file=f)
            for sentinel, anchor_list in self.sentinel_to_anchor.items():
                for anchor in anchor_list:
                    print(
                        f"{sentinel}\t{anchor.snarl_id}\t{anchor.basepairlength}\t{anchor.genomic_position}\t{anchor!r}\t{anchor.bandage_representation()}\t{anchor.get_reference_paths()}",
                        file=f,
                    )
