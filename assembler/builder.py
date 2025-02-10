# bdsg import
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
)
from assembler.rev_c import rev_c
from assembler.node import Node
from assembler.anchor import Anchor

# other imports
import time
# import json
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
        self.keep_path_scan: bool = True
        self.path_orientation: bool = True
        self.snarl_boundaries: list = [dict(), dict()]
        self.current_snarl_start: int = 0
        self.peek_orientations = []
        self.count_in_path: bool = True
        self.snarl_id: int = 0

        self.current_anchor: Anchor = Anchor()
        self.verbose = False
        self.ref_path_name = "CHM13".casefold()
        self.paths_handles = []

        # variables used for debugging
        self.num_usable_bubbles = 0
        self.num_used_bubbles = 0

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

    def get_path_orientation_iteratee(self, step_handle) -> bool:
        node_handle = self.graph.get_handle_of_step(step_handle)
        rev = True
        if self.graph.get_is_reverse(node_handle):
            rev = False
        self.peek_orientations.append(rev)
        if len(self.peek_orientations) >= PEEK_SIZE:
            num_forward = self.peek_orientations.count(True)
            # print(f"Path orientation list: {self.peek_orientations!r}", file=stderr)
            # self.path_orientation = (
            #     FORWARD_DICTIONARY if num_forward >= (PEEK_SIZE // 2) else 1
            # )
            # print(f"Path orientation stored: {self.path_orientation!r}", file=stderr)
            self.peek_orientations = []
            return False
        return True

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
        # print(f"Visiting node {node_id}.", file=stderr)
        if (
            not self.keep_path_scan
            and self.snarl_boundaries[self.path_orientation][self.current_snarl_start][
                0
            ]
            != node_id
        ):
            self.current_anchor.add(
                Node(
                    self.graph.get_id(node_handle),
                    self.graph.get_length(node_handle),
                    not (self.graph.get_is_reverse(node_handle)),
                )
            )
            # print(f"Adding node {self.graph.get_id(node_handle)} to anchor for start {self.current_snarl_start} - {self.snarl_boundaries[self.path_orientation][self.current_snarl_start]}", file=stderr)
            return True

        elif (
            not self.keep_path_scan
            and self.snarl_boundaries[self.path_orientation][self.current_snarl_start][
                0
            ]
            == node_id
        ):
            self.current_anchor.add(
                Node(
                    self.graph.get_id(node_handle),
                    self.graph.get_length(node_handle),
                    not (self.graph.get_is_reverse(node_handle)),
                )
            )
            # print(
            #     f"Found anchor end at {node_id} (start:{self.current_snarl_start}, end:{self.snarl_boundaries[self.path_orientation][self.current_snarl_start][0]})",
            #     file=stderr,
            # )
            #at the end of the snarl. Need to dump the traversal
            # print(
            #     f" Anchor has {len(self.current_anchor)} nodes.",
            #     end=" ",
            #     file=stderr,
            # )
            self.current_anchor.compute_bp_length()
            if (
                self.current_anchor.baseparilength >= MIN_ANCHOR_LENGTH
                and len(self.current_anchor) >= MIN_NODES_IN_ANCHOR
            ):
                sentinel: int = self.current_anchor.get_sentinel_id()
                # print(f"Sentinel is {sentinel}", file=stderr)
                if sentinel not in self.sentinel_to_anchor:
                    # print(f"Added new anchor at sentinel {sentinel}", file=stderr)
                    self.sentinel_to_anchor[sentinel] = [self.current_anchor]

                else:
                    insert = True
                    for inserted_anchor in self.sentinel_to_anchor[sentinel]:
                        # verify that the anchor is not already existing in the dictionary
                        if self.current_anchor == inserted_anchor:
                            insert = False
                    if insert:
                        # print(f"Added +1 anchor at sentinel {sentinel}", file=stderr)
                        self.sentinel_to_anchor[sentinel].append(self.current_anchor)

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
            # print(f"Adding node {self.graph.get_id(node_handle)} to anchor", file=stderr)
            self.current_anchor.add_snarl_id(
                self.snarl_boundaries[self.path_orientation][node_id][1]
            )
            self.keep_path_scan = False
            return True

        # returning True to keep the iteration going
        return True

    def get_edge_snarl(self, snarl_net_handle) -> None:
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
        self.num_usable_bubbles += 1
        start_node_handle, end_node_handle = self.get_snarl_boundaries_handle(
            snarl_net_handle
        )

        snarl_boundaries = (
            (self.graph.get_id(end_node_handle), self.graph.get_id(start_node_handle))
            if self.graph.get_id(end_node_handle) < self.graph.get_id(start_node_handle)
            else (
                self.graph.get_id(start_node_handle),
                self.graph.get_id(end_node_handle),
            )
        )

        #if this node is already used by another bubble, skip it.
        #if self.snarl_boundaries[REVERSE_DICTIONARY].get(snarl_boundaries[0]) != None:
        #    return 
        
        self.snarl_boundaries[FORWARD_DICTIONARY][snarl_boundaries[0]] = (
            snarl_boundaries[1],
            self.snarl_id,
        )
        self.snarl_boundaries[REVERSE_DICTIONARY][snarl_boundaries[1]] = (
            snarl_boundaries[0],
            self.snarl_id,
        )
        self.num_used_bubbles += 1
        return

    def print_anchor_boundaries_dict(self, file_path):
        print(f"Printing to {file_path}.forward_dict.csv")
        with open(f"{file_path}.forward_dict.csv", "w") as f:
            for el in self.snarl_boundaries[FORWARD_DICTIONARY]:
                print(
                    f"{el},{self.snarl_boundaries[FORWARD_DICTIONARY][el][0]}", file=f
                )
        print(f"Printing to {file_path}.reverse_dict.csv")
        with open(f"{file_path}.reverse_dict.csv", "w") as f:
            for el in self.snarl_boundaries[REVERSE_DICTIONARY]:
                print(
                    f"{el},{self.snarl_boundaries[REVERSE_DICTIONARY][el][0]}", file=f
                )

    def collect_path_handles(self, step_handle):
        path_handle = self.graph.get_path_handle_of_step(step_handle)
        self.paths_handles.append(path_handle)  # self.graph.get_path_name()
        return True

    def get_snalrs_from_paths(self) -> None:
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

        #collecting path handles to scan in the graph
        #self.graph.for_each_path_handle(self.collect_path_handles)
        for node in self.snarl_boundaries[0]:
            self.graph.for_each_step_on_handle(self.graph.get_handle(node), self.collect_path_handles )
            # self.graph.for_each_step_on_handle(self.graph.get_handle(end), self.collect_path_handles )

        self.paths_handles = [set(self.paths_handles)]
        
        #scan path handles to obtain the alleles in the snarls.
        print(f"Ready to process {len(self.paths_handles)} paths.")
        for path_handle in self.paths_handles:
            self.current_snarl_start = -1
            self.keep_path_scan = True
            self.count_in_path = True
            self.current_anchor = Anchor()
            print(
                f"Processing path {self.graph.get_path_name(path_handle)}...",
                end=" ",file=stderr,
            )
            t0 = time.time()
            self.graph.for_each_step_in_path(
                path_handle, self.get_path_orientation_iteratee
            )
 
            self.graph.for_each_step_in_path(path_handle, self.traverse_step_iteratee)
            print(f"in {time.time()-t0:.2f} seconds.")

    def generate_anchors_boundaries(self):
        for _, snarl_net_handle in enumerate(self.leaf_snarls):
            self.get_edge_snarl(snarl_net_handle)
        print(f"DICTIONARY",file=stderr)
        for key, value in self.snarl_boundaries[FORWARD_DICTIONARY].items():
            print(f"{key}\t{value}", file=stderr)
        # print(
        #     f"# leaf snarls: {len(self.leaf_snarls)} | # element in start snarl boundaries: {len(self.snarl_boundaries[FORWARD_DICTIONARY])}, end: {len(self.snarl_boundaries[REVERSE_DICTIONARY])}"
        # )

    def fill_anchor_dictionary(self) -> None:
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
        self.generate_anchors_boundaries()
        print(
            f"Snarl Boundaries computed in {time.time()-t1:.2f}",
            file=stderr,
        )

        # for idx, snarl_net_handle in enumerate(self.leaf_snarls):
        t2 = time.time()
        self.get_snalrs_from_paths()
        # self.fill_anchor_dictionary_single_snarl(snarl_net_h)
        print(
            f"Snarl dictionary computed in {time.time()-t2:.2f}. Total time: {time.time()-t0:.2f}.",
            file=stderr,
        )
        print(f"Num used bubbles: {self.num_used_bubbles} ; Num usable bubbles: {self.num_usable_bubbles} ; ratio {self.num_used_bubbles/self.num_usable_bubbles}.",file=stderr)

    ### HELPER FUNCTIONS ###

    def get_nodes_in_snarl(self, rnarl_net_handle) -> list:
        nodes_inside = []
        self.index.traverse_decomposition_helper(
        rnarl_net_handle,  
        snarl_iteratee=lambda s: True,
        chain_iteratee=lambda c: True,
        node_iteratee=lambda n: (self.graph.get_id(nodes_inside.append(self.index.get_handle(n, self.graph))) or True)
        )


    def get_left_neighbours(self, node_handle) -> list :
        left_neighbors = []
        self.graph.follow_edges(node_handle, go_left=True, iteratee=lambda neighbor: left_neighbors.append(neighbor) or True)
        return left_neighbors
    
    def get_right_neighbours(self, node_handle) -> list :
        right_neighbors = []
        self.graph.follow_edges(node_handle, go_left=False, iteratee=lambda neighbor: right_neighbors.append(neighbor) or True)
        return right_neighbors

    # Returns a tuple containing the start and end boundaries of a snarl
    def get_snarl_boundaries_handle_chuncked_graph(self, snarl_net_handle) -> tuple:
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

        nodes_in_snarl = self.get_nodes_in_snarl(snarl_net_handle)

        start_node_handle = self.index.get_handle(start_bound, self.graph)
        end_node_handle = self.index.get_handle(end_bound, self.graph)
        
        start_node_neighbours = self.get_left_neighbours(start_node_handle)
        end_node_neighbours = self.get_right_neighbours(end_node_handle)

        set_next_start_node = False
        if len(start_node_neighbours) == 1:
            start_node_neighbour = self.graph.get_id(start_node_neighbours[0])
            if start_node_neighbour not in nodes_in_snarl:
                start_node_handle = start_node_neighbours[0]
                set_next_start_node = True
        if not set_next_start_node:
            start_node_neighbours = self.get_left_neighbours(start_node_handle)
            if len(start_node_neighbours) == 1:
                start_node_neighbour = self.graph.get_id(start_node_neighbours[0])
                if start_node_neighbour not in nodes_in_snarl:
                    start_node_handle = start_node_neighbours[0]

        if len(end_node_neighbours) == 1:
            end_node_handle = end_node_neighbours[0]
        

        boundary = (
            start_node_handle,
            end_node_handle,
        )
        return boundary

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
        path_names = self.get_path_names()

        for path in path_names:
            if self.ref_path_name in path.casefold():
                graph_path_name = path

        if self.graph.has_path(graph_path_name):
            print(f"Found path {graph_path_name}", file=stderr)
            path_handle = self.graph.get_path_handle(graph_path_name)
        else:
            print(f"WARNING: Could not find CHM13 path in graph", file=stderr)
            return

        # generate main path (node_id,lenght_in_chm13)
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
            print(f"Sentinel_node\tAnchor_length\tAnchor_pos_in_ref_path\tAnchor_path\tAnchor_nodes_copypaste_bandage",file=f)
            for sentinel, anchor_list in self.sentinel_to_anchor.items():
                for anchor in anchor_list:
                    print(
                        f"{sentinel}\t{anchor.baseparilength}\t{anchor.genomic_position}\t{anchor!r}\t{anchor.bandage_representation()}",
                        file=f,
                    )
    
    def print_paths_used(self, out_f) -> None:
        with open(out_f, "w") as f:
            for path in self.paths_handles:
                print(f"{self.graph.get_path_name(path)}", file=f)