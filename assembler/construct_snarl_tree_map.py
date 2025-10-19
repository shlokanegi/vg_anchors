"""
Constructs a hierarchical snarl tree map from a variation graph.

This script traverses a snarl decomposition of a variation graph, calculating
the number of leaf snarls in the subtree of each snarl and chain. It supports
parallel processing for large graphs by dividing the traversal work at a
predefined depth in the tree.

The output is a JSON file representing the snarl tree, where each key is a
snarl or chain ID and the value contains its children and the total count
of leaf snarls in its subtree.
"""
import os
import multiprocessing
from bdsg.bdsg import PackedGraph
from bdsg.bdsg import SnarlDistanceIndex
import time
import json
from sys import stderr
import argparse

PARALLELIZATION_DEPTH = 2

g_snarl_decomposition_obj = None

def init_worker(snarl_decomposition_obj):
    """
    Initializes each worker process for parallel execution.

    This function is called by multiprocessing.Pool to set up the global
    SnarlDecomposition object in each worker, allowing for memory sharing
    of the loaded graph and index via copy-on-write.
    """
    global g_snarl_decomposition_obj
    g_snarl_decomposition_obj = snarl_decomposition_obj

def process_subtree_worker(args):
    """
    Performs snarl tree traversal on a specific subtree.

    This is the main worker function called by the multiprocessing pool. It
    uses the globally shared SnarlDecomposition object to process a single
    subtree and returns the resulting map.
    """
    worker_id, subtree_root_handle = args
    start_time = time.time()
        
    result = g_snarl_decomposition_obj.traverse_decomposition(subtree_root_handle)
        
    end_time = time.time()
    print(f"[Worker {worker_id}] Finished in {end_time - start_time:.2f} seconds.", file=stderr)
        
    return result

class SnarlDecomposition:
    """
    Manages the decomposition of a variation graph into a snarl tree.

    This class encapsulates the graph and its distance index, providing methods
    to traverse the snarl decomposition and build a hierarchical map.
    """

    def __init__(self):
        """Initializes the SnarlDecomposition object."""
        self.graph = PackedGraph()
        self.index = SnarlDistanceIndex()
        self.snarl_tree_map = {}   # snarl_tree_map[snarl_id] = {[children_ids], count_of_leaf_snarls_in_subtree}

    def build(self, graph_path, index_path):
        """Loads the graph and snarl distance index from disk."""
        t0 = time.time()
        self.graph.deserialize(graph_path)
        print(f"Graph deserialized in {time.time()-t0:.2f}", file=stderr, flush=True)
        self.index.deserialize(index_path)
        print(f"Index deserialized in {time.time()-t0:.2f}", file=stderr, flush=True)


    def get_id(self, net):
        """Generates a unique, human-readable ID for a net object."""
        if self.index.is_chain(net) or self.index.is_snarl(net):
            left_bound_net = self.index.get_bound(net, False, False)
            right_bound_net = self.index.get_bound(net, True, False)
            left_bound_node_id = self.graph.get_id(self.index.get_handle(left_bound_net, self.graph))
            right_bound_node_id = self.graph.get_id(self.index.get_handle(right_bound_net, self.graph))
            if left_bound_node_id < right_bound_node_id:
                left_bound_node_id, right_bound_node_id = right_bound_node_id, left_bound_node_id
            bound_str = str(left_bound_node_id) + "-" + str(right_bound_node_id)
            if self.index.is_chain(net):
                return "c" + bound_str
            else:
                return "s" + bound_str
        elif self.index.is_root(net):
            return "root"
        elif self.index.is_node(net):
            return str(self.graph.get_id(self.index.get_handle(net, self.graph)))
        else:
            raise ValueError(f"Net {net} is not a chain, snarl, node or root")

    def check_snarl_in_children_iteratee(self, child_net_handle) -> bool:
        """Callback function to check if a net has snarls as children."""
        if self.index.is_snarl(child_net_handle):
            self.contains_child_snarls = True
            return False
        return True

    def is_leaf_snarl(self, net):
        """Determines if a snarl is a leaf in the snarl tree."""
        self.contains_child_snarls = False

        snarl_children: list = []
        self.index.for_each_child(
            net, lambda y: snarl_children.append(y) or True
        )

        for snarl_child in snarl_children:
            self.index.for_each_child(snarl_child, self.check_snarl_in_children_iteratee)

        if not self.contains_child_snarls:
            return True

    
    def traverse_decomposition_helper(self, handle, snarl_tree_map):
        """
        Recursively traverses a subtree to count leaf snarls.

        This is the core recursive function for building the snarl tree map.
        It performs a depth-first traversal, counting leaf snarls and
        propagating the counts up the tree.
        """
        map_key = self.get_id(handle)
        
        # The stopping condition is when the current net is a leaf snarl
        if self.index.is_snarl(handle) and self.is_leaf_snarl(handle):
            snarl_tree_map[map_key] = ([], 1)
            return 1, map_key
        elif self.index.is_node(handle):
            return 0, map_key

        cnt_leaf_snarls_in_subtree = 0
        children_handles = []
        children_map_keys = []
        self.index.for_each_child(handle, lambda y: (children_handles.append(y) if not self.index.is_node(y) else True) or True)
        for child_handle in children_handles:
            cnt_leaf_snarls_in_child_subtree, child_map_key = self.traverse_decomposition_helper(child_handle, snarl_tree_map)
            children_map_keys.append(child_map_key)
            # Once the recursive call for a child returns, it will have returned the total number of leaf snarls in that child's entire subtree. This count is added to the current snarl's running total.
            cnt_leaf_snarls_in_subtree += cnt_leaf_snarls_in_child_subtree
        
        # After the for loop ends, the function has visited all the current snarl's children and has received the total leaf snarl counts for each of them.
        snarl_tree_map[map_key] = (children_map_keys, cnt_leaf_snarls_in_subtree)
        
        return cnt_leaf_snarls_in_subtree, map_key


    def traverse_decomposition(self, subtree_root_handle):
        """
        Traverses a subtree to build a snarl tree map.

        This function is re-entrant and can be called from multiple processes.
        It initializes a new map for each traversal and returns it.
        """
        snarl_tree_map = {}
        self.traverse_decomposition_helper(subtree_root_handle, snarl_tree_map)
        # Handle the case where the root of the subtree is itself a leaf snarl
        if not snarl_tree_map and self.index.is_snarl(subtree_root_handle) and self.is_leaf_snarl(subtree_root_handle):
             map_key = self.get_id(subtree_root_handle)
             snarl_tree_map[map_key] = ([], 1)
        return snarl_tree_map


    def collect_handles_at_depth(self, depth):
        """
        Collects all snarl and chain handles at a specific depth.

        This function now works by iterating through the nodes of parent
        handles to find the unique child snarls and chains.
        """
        if depth == 0:
            return [self.index.get_root()]

        parent_handles = self.collect_handles_at_depth(depth - 1)

        child_handles_at_parallelization_depth = []

        for parent_handle in parent_handles:
            self.index.for_each_child(parent_handle, lambda y: (child_handles_at_parallelization_depth.append(y) if not self.index.is_node(y) else True) or True)

        return child_handles_at_parallelization_depth


    def build_top_tree(self, handle, current_depth):
        """
        Recursively reconstructs the top levels of the snarl tree map.

        After parallel processing of the lower levels, this function is called
        on the root to rebuild the top of the tree, integrating the results
        from the worker processes.
        """
        if current_depth >= PARALLELIZATION_DEPTH:
            return self.snarl_tree_map.get(self.get_id(handle), ([], 0))[1]

        map_key = self.get_id(handle)
        
        if map_key in self.snarl_tree_map and self.snarl_tree_map[map_key][1] > 0:
            return self.snarl_tree_map[map_key][1]

        children_handles = []
        self.index.for_each_child(handle, lambda y: (children_handles.append(y) if not self.index.is_node(y) else True) or True)

        children_map_keys = []
        total_leaf_snarls = 0

        for child_handle in children_handles:
            children_map_keys.append(self.get_id(child_handle))
            count = self.build_top_tree(child_handle, current_depth + 1)
            total_leaf_snarls += count

        if not children_handles and self.index.is_snarl(handle) and self.is_leaf_snarl(handle):
            total_leaf_snarls = 1

        self.snarl_tree_map[map_key] = (children_map_keys, total_leaf_snarls)
        return total_leaf_snarls


def main():
    """
    Main function to orchestrate the snarl tree map construction.

    Parses command-line arguments, loads the graph and index, manages the
    parallel processing of subtrees, reconstructs the final tree, and saves
    the result to a JSON file.
    """
    parser = argparse.ArgumentParser(description="Construct a snarl tree map from a variation graph.")
    parser.add_argument("-g", "--graph", required=True, help="Path to the variation graph (.pg file).")
    parser.add_argument("-i", "--index", required=True, help="Path to the snarl distance index (.dist file).")
    parser.add_argument("-o", "--output-json", default="snarl_tree_map.json", help="Path to the output JSON file.")
    args = parser.parse_args()
    
    print("Starting parallel snarl tree map construction.", file=stderr)

    main_obj = SnarlDecomposition()
    main_obj.build(args.graph, args.index)

    handles_to_process = main_obj.collect_handles_at_depth(PARALLELIZATION_DEPTH)
    print(f"Found {len(handles_to_process)} handles to process.", file=stderr)
    
    num_processes = min(multiprocessing.cpu_count(), len(handles_to_process)) if handles_to_process else 0
    
    print(f"Running with {num_processes} parallel processes.", file=stderr)

    if num_processes > 0:
        # We use the 'fork' start method (default on Linux) to allow for memory
        # sharing of the main_obj between processes (copy-on-write).
        with multiprocessing.Pool(processes=num_processes, initializer=init_worker, initargs=(main_obj,)) as pool:
            worker_args = [(i, handle) for i, handle in enumerate(handles_to_process)]
            results = pool.map(process_subtree_worker, worker_args)
            
            merged_map = {}
            for result_map in results:
                merged_map.update(result_map)
            
            # The main object's map now contains the results from all workers
            main_obj.snarl_tree_map = merged_map
    else: # If no handles to process, the tree is smaller than PARALLELIZATION_DEPTH
        # If the tree is very small, just traverse it sequentially.
        main_obj.snarl_tree_map = main_obj.traverse_decomposition(main_obj.index.get_root())

    # Reconstruct the top levels of the tree using the processed subtrees
    main_obj.build_top_tree(main_obj.index.get_root(), 0)
    
    # save to json
    with open(args.output_json, "w") as f:
        json.dump(main_obj.snarl_tree_map, f, indent=4)

if __name__ == "__main__":
    main()
