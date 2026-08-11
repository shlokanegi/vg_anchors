# bdsg import
# import json
import time
from bdsg.bdsg import SnarlDistanceIndex
from bdsg.bdsg import PackedGraph

# package import
from assembler.config import settings
from assembler.node import Node
from assembler.anchor import Anchor
from assembler.bdsg_compat import check_libbdsg_index_support, explain_traversal_failure

# other imports
import time
from sys import stderr
import pickle


class ChunkerAnchorDictionary:
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
        self.step_counts = dict()

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
        check_libbdsg_index_support()
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
        try:
            self.index.traverse_decomposition(
                self.check_leaf_snarl_iteratee,  # snarl_iteratee
                lambda x: True,  #  chain_iteratee
                lambda y: True,  # node_iteratee
            )
        except RuntimeError as error:
            explanation = explain_traversal_failure(error)
            if explanation is None:
                raise
            raise RuntimeError(explanation) from error
        return None


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

        if settings.DEBUG:
            print(f"In current path, traversing node id {node_id}", end="\n")

        if (
            not self.keep_path_scan
            and self.snarl_boundaries[self.path_orientation][self.current_snarl_start][settings.END_NODE_POS] != node_id
            and node_id in self.snarl_boundaries[self.path_orientation][self.current_snarl_start][2]
        ):
            self.current_anchor.add(
                Node(
                    self.graph.get_id(node_handle),
                    self.graph.get_length(node_handle),
                    not (self.graph.get_is_reverse(node_handle)),
                )
            )
            if settings.DEBUG:
                print(f"Adding node {self.graph.get_id(node_handle)} to anchor gets {self.current_anchor!r}", end="\n")

            return True

        elif (
            not self.keep_path_scan
        ):
            if self.snarl_boundaries[self.path_orientation][self.current_snarl_start][settings.END_NODE_POS] == node_id:
                self.current_anchor.add(
                    Node(
                        self.graph.get_id(node_handle),
                        self.graph.get_length(node_handle),
                        not (self.graph.get_is_reverse(node_handle)),
                    )
                )
                
                self.current_anchor.compute_snarl_boundary()                
                self.current_anchor.path_orientation = self.current_anchor._nodes[0].id < self.current_anchor._nodes[-1].id
                
                if (
                    len(self.current_anchor) >= settings.MIN_NODES_IN_ANCHOR
                ):
                    sentinel: int = self.current_anchor.get_sentinel_id()
                    if sentinel not in self.sentinel_to_anchor:
                        self.current_anchor.add_reference_path(self.curr_path_name)
                        self.sentinel_to_anchor[sentinel] = [self.current_anchor]
                        if settings.DEBUG:
                            print(f"Final anchor is {self.current_anchor!r} whose sentinal is {sentinel} and length {self.current_anchor.basepairlength}")

                    else:
                        insert = True
                        for id, inserted_anchor in enumerate(self.sentinel_to_anchor[sentinel]):
                            # verify that the anchor is not already existing in the dictionary
                            if self.current_anchor == inserted_anchor:  # if current anchor is same as the one already in the list at index 'id', then just update the path variable of the anchor
                                self.sentinel_to_anchor[sentinel][id].add_reference_path(self.curr_path_name)
                                if settings.DEBUG:
                                    print(f"Final anchor is {self.current_anchor!r} whose sentinal is {sentinel} and length {self.current_anchor.basepairlength}")
                                insert = False
                        if insert:
                            # but, if current anchor is not already in the list, then add it to the list and also update path variable
                            self.current_anchor.add_reference_path(self.curr_path_name)
                            self.sentinel_to_anchor[sentinel].append(self.current_anchor)
                            if settings.DEBUG:
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
            if settings.DEBUG:
                print(f"Adding node {self.graph.get_id(node_handle)} to anchor. Corresponding boundary node is {self.snarl_boundaries[self.path_orientation][self.current_snarl_start][settings.END_NODE_POS]}", end="\n")
            self.current_anchor.add_snarl_id(
                self.snarl_boundaries[self.path_orientation][node_id][settings.SNARL_ID_POS]
            )
            self.used_bubbles[self.snarl_boundaries[self.path_orientation][node_id][settings.SNARL_ID_POS]] = True
            self.keep_path_scan = False
            return True

        # returning True to keep the iteration going
        return True


    def get_edge_snarl(self, snarl_net_handle) -> None:
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
        
        start_node_handle, end_node_handle, nodes_to_select = self.get_snarl_boundaries_handle(snarl_net_handle)

        snarl_boundary = (self.graph.get_id(end_node_handle), self.graph.get_id(start_node_handle))

        # # select one node per allele in the leaf snarl
        # nodes_to_select = []
        # for chain_id, nodes in chain_to_sentinel_nodes_dict.items():
        #     nodes_to_select.append(nodes[0])
        
        self.snarl_boundaries[settings.FORWARD_DICTIONARY][snarl_boundary[0]] = (
            snarl_boundary[1],
            self.num_usable_bubbles,    # SNARL ID
            nodes_to_select
        )

        return

    def print_anchor_boundaries_dict(self, file_path):
        if settings.DEBUG:
            print(f"Printing to {file_path}.forward_dict.csv")
        with open(f"{file_path}.forward_dict.csv", "w") as f:
            for el in self.snarl_boundaries[settings.FORWARD_DICTIONARY]:
                print(
                    f"{el},{self.snarl_boundaries[settings.FORWARD_DICTIONARY][el][settings.END_NODE_POS]},{self.snarl_boundaries[settings.FORWARD_DICTIONARY][el][2]}", file=f
                )
        if settings.DEBUG:
            print(f"Printing to {file_path}.reverse_dict.csv")
        with open(f"{file_path}.reverse_dict.csv", "w") as f:
            for el in self.snarl_boundaries[settings.REVERSE_DICTIONARY]:
                print(
                    f"{el},{self.snarl_boundaries[settings.REVERSE_DICTIONARY][el][settings.END_NODE_POS]},{self.snarl_boundaries[settings.REVERSE_DICTIONARY][el][2]}", file=f
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

        if settings.DEBUG:
            print(f"TOT PATHS COLLECTED: {len(self.path_names)}")
        self.path_names = sorted(list(set(self.path_names)))
        #scan path handles to obtain the alleles in the snarls.
        if settings.DEBUG:
            print(f"Ready to process {len(self.path_names)} paths...", end = ' ')
        t_0 = time.time()
        for path_name in self.path_names:
            if settings.DEBUG:
                print(f"Processing path {path_name}...")
            for path_orientation in [settings.REVERSE_DICTIONARY, settings.FORWARD_DICTIONARY]:
                print(f"Processing path {path_name} in orientation {path_orientation}")
                path_handle = self.graph.get_path_handle(path_name)
                self.curr_path_name = path_name
                if settings.DEBUG:
                    print(f"Currently processing path {self.curr_path_name}...", end="\n")
                self.current_snarl_start = -1
                self.keep_path_scan = True
                self.count_in_path = True
                self.current_anchor = Anchor()
                self.peek_orientations = []
                self.path_orientation=path_orientation

                if settings.DEBUG:
                    print(f"With path_orientation {self.path_orientation}", end="\n")
                self.graph.for_each_step_in_path(path_handle, self.traverse_step_iteratee)

            if settings.DEBUG:
                print(f"done in {time.time()-t_0}")


    def generate_anchors_boundaries(self):
        """
        This function sorts leaf snarl handle list based on snarl orientation, so that all snarl handles are in ascending order of occurrence.
        It also updates the snarl dictionaries
        """
        
        for _, snarl_net_handle in enumerate(self.leaf_snarls):
            self.get_edge_snarl(snarl_net_handle)


    def get_step_counts_from_sentinel_nodes_of_snarls(self) -> None:
        """
        This function gets the step counts from the sentinel nodes of the snarls.
        """
        chk = 0
        t_0 = time.time()
        for start_bound, (end_bound, snarl_id, nodes_inside) in self.snarl_boundaries[settings.FORWARD_DICTIONARY].items():
            for node_id in nodes_inside:
                step_counts = 0
                
                def step_count_callable(step_handle):
                    nonlocal step_counts
                    step_counts += 1
                    return True
                
                node_handle = self.graph.get_handle(node_id)
                self.graph.for_each_step_on_handle(node_handle, step_count_callable)
                
                # Update step count dictionary
                snarl_bound_id = f"{start_bound}_{end_bound}"
                if snarl_bound_id not in self.step_counts:
                    self.step_counts[snarl_bound_id] = []
                self.step_counts[snarl_bound_id].append(step_counts)
                self.flag_snarl_bound_id = snarl_bound_id
            
            chk += 1
            if chk % 100000 == 0:
                # print the time taken to process the snarls
                print(f"Time taken to process {chk} snarls: {time.time() - t_0:.2f} seconds", flush=True, file=stderr)
                print(f"Last snarl ID inserted is {self.flag_snarl_bound_id} and the number of steps is {self.step_counts[self.flag_snarl_bound_id]}")


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
        self.generate_anchors_boundaries()
        print(
            f"Snarl Boundaries computed in {time.time()-t1:.2f}",
            file=stderr,
        )

        t2 = time.time()
        self.get_step_counts_from_sentinel_nodes_of_snarls()

        print(
            f"Snarl dictionary computed in {time.time()-t2:.2f}. Total time: {time.time()-t0:.2f}.",
            file=stderr,
        )
        num_used_bubbles = len([_ for _ in self.used_bubbles if self.used_bubbles[_] == 1])
        print(f"Num used bubbles: {num_used_bubbles} ; Num usable bubbles: {self.num_usable_bubbles} ; ratio {num_used_bubbles/self.num_usable_bubbles}.",file=stderr)


    ### HELPER FUNCTIONS ###
    def get_sentinel_nodes_in_snarl(self, snarl_net_handle) -> list:

        nodes_to_select = []
        
        def chain_iteratee_callable(n):
            # make chain ID
            left_bound_net = self.index.get_bound(n, False, False)
            left_bound_node_id = self.graph.get_id(self.index.get_handle(left_bound_net, self.graph))
            nodes_to_select.append(left_bound_node_id)
            return True

        self.index.for_each_child(snarl_net_handle, chain_iteratee_callable)

        return nodes_to_select

        # chain_to_nodes_dict = {}

        # def node_iteratee_callable(n):
        #     parent_chain_handle = self.index.get_parent(n)
        #     # make chain ID
        #     left_bound_net = self.index.get_bound(parent_chain_handle, False, False)
        #     right_bound_net = self.index.get_bound(parent_chain_handle, True, False)
        #     left_bound_node_id = self.graph.get_id(self.index.get_handle(left_bound_net, self.graph))
        #     right_bound_node_id = self.graph.get_id(self.index.get_handle(right_bound_net, self.graph))
        #     if left_bound_node_id < right_bound_node_id:
        #         left_bound_node_id, right_bound_node_id = right_bound_node_id, left_bound_node_id
        #     bound_str = str(left_bound_node_id) + "-" + str(right_bound_node_id)
        #     chain_id = "c" + bound_str

        #     if chain_id not in chain_to_nodes_dict:
        #         chain_to_nodes_dict[chain_id] = []
        #     chain_to_nodes_dict[chain_id].append((self.graph.get_id(self.index.get_handle(n, self.graph))))

        #     return True

        # self.index.traverse_decomposition_helper(
        #     snarl_net_handle,  
        #     snarl_iteratee=lambda s: True,  # Ignore snarls
        #     chain_iteratee=lambda c: True,  # Ignore chains
        #     node_iteratee=node_iteratee_callable
        # )

        # return chain_to_nodes_dict  # Return the collected node IDs

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
        nodes_to_select = self.get_sentinel_nodes_in_snarl(snarl_net_handle)

        return (start_bound_handle,end_bound_handle,nodes_to_select)


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

    def print_hap_counts(self, out_f) -> None:
        with open(out_f, "w") as f:
            print(f"snarl_id\tstart_bound\tend_bound\tstep_counts",file=f)

            for snarl_bound_id, step_counts_list in self.step_counts.items():
                start_bound, end_bound = snarl_bound_id.split("_")
                print(f"{snarl_bound_id}\t{start_bound}\t{end_bound}\t{','.join(map(str, step_counts_list))}", file=f)
