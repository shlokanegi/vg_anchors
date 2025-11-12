/**
 * Constructs a hierarchical snarl tree map from a variation graph.
 * 
 * This program traverses a snarl decomposition of a variation graph, calculating
 * the number of leaf snarls in the subtree of each snarl and chain. It supports
 * parallel processing for large graphs using C++ multithreading.
 * 
 * The output is a JSON file representing the snarl tree, where each key is a
 * snarl or chain ID and the value contains its children and the total count
 * of leaf snarls in its subtree.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <mutex>
#include <chrono>
#include <algorithm>
#include <memory>
#include <sstream>
#include <iomanip>

#include "bdsg/packed_graph.hpp"
#include "bdsg/snarl_distance_index.hpp"
#include "handlegraph/types.hpp"
#include "handlegraph/util.hpp"

using namespace std;
using namespace bdsg;
using namespace handlegraph;

const int PARALLELIZATION_DEPTH = 2;

// Structure to hold the tree map entry: (children_ids, leaf_count)
struct TreeMapEntry {
    vector<string> children_ids;
    int leaf_count;
};

// Structure to hold chunk point: (start_node_id, end_node_id, leaf_snarl_count, path_name, orientations)
struct ChunkPoint {
    nid_t start_node;
    nid_t end_node;
    int leaf_count;
    string path_name;
    bool start_is_reverse;  // true if start node is reverse orientation
    bool end_is_reverse;    // true if end node is reverse orientation
    string root_chain_id;   // the ID of the root chain
};

class SnarlTreeBuilder {
private:
    PackedGraph graph;
    SnarlDistanceIndex index;
    map<string, TreeMapEntry> snarl_tree_map;
    vector<ChunkPoint> chunk_points;  // List of chunk points for parallel processing
    mutex map_mutex;  // For thread-safe map operations
    int target_leaf_snarls_per_chunk;  // Target number of leaf snarls per chunk
    
    /**
     * Gets the boundary node IDs of a chain or snarl.
     * Returns a pair of (start_node_id, end_node_id).
     */
    pair<nid_t, nid_t> get_boundary_nodes(const net_handle_t& net) {
        net_handle_t left_bound_net = index.get_bound(net, false, false);
        net_handle_t right_bound_net = index.get_bound(net, true, false);
        nid_t left_bound_node_id = graph.get_id(index.get_handle(left_bound_net, &graph));
        nid_t right_bound_node_id = graph.get_id(index.get_handle(right_bound_net, &graph));
        
        // Return in the order they appear (left, right)
        return {left_bound_node_id, right_bound_node_id};
    }
    
    /**
     * Gets the boundary node IDs and orientations of a chain or snarl.
     * Returns a tuple of (start_node_id, end_node_id, start_is_reverse, end_is_reverse).
     */
    tuple<nid_t, nid_t, bool, bool> get_boundary_nodes_with_orientation(const net_handle_t& net) {
        net_handle_t left_bound_net = index.get_bound(net, false, false);
        net_handle_t right_bound_net = index.get_bound(net, true, false);
        
        handle_t left_handle = index.get_handle(left_bound_net, &graph);
        handle_t right_handle = index.get_handle(right_bound_net, &graph);
        
        nid_t left_bound_node_id = graph.get_id(left_handle);
        nid_t right_bound_node_id = graph.get_id(right_handle);
        
        bool left_is_reverse = graph.get_is_reverse(left_handle);
        bool right_is_reverse = graph.get_is_reverse(right_handle);
        
        // Return in the order they appear (left, right)
        return {left_bound_node_id, right_bound_node_id, left_is_reverse, right_is_reverse};
    }
    
    /**
     * Gets the path name (CHM13 reference contig) that contains the given node.
     * Returns the first path found that traverses this node.
     * If no path is found, returns "unknown".
     */
    string get_path_for_node(nid_t node_id) {
        string path_name = "unknown";
        
        // Get the handle for this node
        handle_t node_handle = graph.get_handle(node_id);
        
        // Iterate through all paths that traverse this node
        graph.for_each_step_on_handle(node_handle, [&](const step_handle_t& step) {
            // Get the path handle from this step
            path_handle_t path = graph.get_path_handle_of_step(step);
            
            // Get the path name
            string current_path_name = graph.get_path_name(path);
            
            // Prefer paths that contain "CHM13" or "GRCh38" or similar reference indicators
            // If this is a reference path, use it
            if (current_path_name.find("CHM13") != string::npos || 
                current_path_name.find("GRCh38") != string::npos ||
                current_path_name.find("chm13") != string::npos ||
                current_path_name.find("grch38") != string::npos) {
                path_name = current_path_name;
                return false;  // Stop iteration once we find a reference path
            }
            
            // Otherwise, if we haven't found anything yet, use this path
            if (path_name == "unknown") {
                path_name = current_path_name;
            }
            
            return true;  // Continue iteration
        });
        
        return path_name;
    }
    
    /**
     * Generates a unique, human-readable ID for a net object.
     */
    string get_id(const net_handle_t& net) {
        if (index.is_chain(net) || index.is_snarl(net)) {
            net_handle_t left_bound_net = index.get_bound(net, false, false);
            net_handle_t right_bound_net = index.get_bound(net, true, false);
            nid_t left_bound_node_id = graph.get_id(index.get_handle(left_bound_net, &graph));
            nid_t right_bound_node_id = graph.get_id(index.get_handle(right_bound_net, &graph));
            
            if (left_bound_node_id < right_bound_node_id) {
                swap(left_bound_node_id, right_bound_node_id);
            }
            
            string bound_str = to_string(left_bound_node_id) + "-" + to_string(right_bound_node_id);
            
            if (index.is_chain(net)) {
                return "c" + bound_str;
            } else {
                return "s" + bound_str;
            }
        } else if (index.is_root(net)) {
            return "root";
        } else if (index.is_node(net)) {
            return to_string(graph.get_id(index.get_handle(net, &graph)));
        } else {
            throw runtime_error("Net is not a chain, snarl, node or root");
        }
    }
    
    /**
     * Determines if a snarl is a leaf in the snarl tree.
     */
    bool is_leaf_snarl(const net_handle_t& net) {
        bool contains_child_snarls = false;
        
        vector<net_handle_t> snarl_children;
        index.for_each_child(net, [&](const net_handle_t& child) {
            snarl_children.push_back(child);
        });
        
        for (const auto& snarl_child : snarl_children) {
            index.for_each_child(snarl_child, [&](const net_handle_t& grandchild) {
                if (index.is_snarl(grandchild)) {
                    contains_child_snarls = true;
                    return false;  // Stop iteration
                }
                return true;
            });
            if (contains_child_snarls) {
                break;
            }
        }
        
        return !contains_child_snarls;
    }
    
    /**
     * Recursively traverses a subtree to count leaf snarls.
     */
    pair<int, string> traverse_decomposition_helper(const net_handle_t& handle, 
                                                     map<string, TreeMapEntry>& local_map) {
        string map_key = get_id(handle);
        
        // Stopping condition: leaf snarl
        if (index.is_snarl(handle) && is_leaf_snarl(handle)) {
            local_map[map_key] = {{}, 1};
            return {1, map_key};
        } else if (index.is_node(handle)) {
            return {0, map_key};
        }
        
        int cnt_leaf_snarls_in_subtree = 0;
        vector<net_handle_t> children_handles;
        vector<string> children_map_keys;
        
        // Collect non-node children
        index.for_each_child(handle, [&](const net_handle_t& child) {
            if (!index.is_node(child)) {
                children_handles.push_back(child);
            }
        });
        
        for (const auto& child_handle : children_handles) {
            auto [cnt_leaf_snarls_in_child, child_map_key] = 
                traverse_decomposition_helper(child_handle, local_map);
            children_map_keys.push_back(child_map_key);
            cnt_leaf_snarls_in_subtree += cnt_leaf_snarls_in_child;
        }
        
        local_map[map_key] = {children_map_keys, cnt_leaf_snarls_in_subtree};
        
        return {cnt_leaf_snarls_in_subtree, map_key};
    }
    
    /**
     * Traverses a subtree to build a snarl tree map.
     * Thread-safe for parallel execution.
     */
    map<string, TreeMapEntry> traverse_decomposition(const net_handle_t& subtree_root_handle) {
        map<string, TreeMapEntry> local_map;
        traverse_decomposition_helper(subtree_root_handle, local_map);
        
        // Handle the case where the root is itself a leaf snarl
        if (local_map.empty() && index.is_snarl(subtree_root_handle) && 
            is_leaf_snarl(subtree_root_handle)) {
            string map_key = get_id(subtree_root_handle);
            local_map[map_key] = {{}, 1};
        }
        
        return local_map;
    }
    
    /**
     * Collects all snarl and chain handles at a specific depth.
     */
    vector<net_handle_t> collect_handles_at_depth(int depth) {
        if (depth == 0) {
            return {index.get_root()};
        }
        
        vector<net_handle_t> parent_handles = collect_handles_at_depth(depth - 1);
        vector<net_handle_t> child_handles;
        
        for (const auto& parent_handle : parent_handles) {
            index.for_each_child(parent_handle, [&](const net_handle_t& child) {
                if (!index.is_node(child)) {
                    child_handles.push_back(child);
                }
            });
        }
        
        return child_handles;
    }
    
    /**
     * Worker function for parallel processing.
     */
    void process_subtree_worker(int worker_id, const net_handle_t& handle) {
        auto local_map = traverse_decomposition(handle);
        
        // Merge into shared map (thread-safe)
        lock_guard<mutex> lock(map_mutex);
        snarl_tree_map.insert(local_map.begin(), local_map.end());
    }
    
    /**
     * Recursively reconstructs the top levels of the snarl tree map.
     */
    int build_top_tree(const net_handle_t& handle, int current_depth) {
        if (current_depth >= PARALLELIZATION_DEPTH) {
            string key = get_id(handle);
            auto it = snarl_tree_map.find(key);
            if (it != snarl_tree_map.end()) {
                return it->second.leaf_count;
            }
            return 0;
        }
        
        string map_key = get_id(handle);
        
        // Check if already computed
        auto it = snarl_tree_map.find(map_key);
        if (it != snarl_tree_map.end() && it->second.leaf_count > 0) {
            return it->second.leaf_count;
        }
        
        vector<net_handle_t> children_handles;
        index.for_each_child(handle, [&](const net_handle_t& child) {
            if (!index.is_node(child)) {
                children_handles.push_back(child);
            }
        });
        
        vector<string> children_map_keys;
        int total_leaf_snarls = 0;
        
        for (const auto& child_handle : children_handles) {
            children_map_keys.push_back(get_id(child_handle));
            int count = build_top_tree(child_handle, current_depth + 1);
            total_leaf_snarls += count;
        }
        
        if (children_handles.empty() && index.is_snarl(handle) && is_leaf_snarl(handle)) {
            total_leaf_snarls = 1;
        }
        
        snarl_tree_map[map_key] = {children_map_keys, total_leaf_snarls};
        return total_leaf_snarls;
    }
    
    /**
     * Escapes a string for JSON output.
     */
    string json_escape(const string& s) {
        ostringstream o;
        for (auto c : s) {
            switch (c) {
                case '"': o << "\\\""; break;
                case '\\': o << "\\\\"; break;
                case '\b': o << "\\b"; break;
                case '\f': o << "\\f"; break;
                case '\n': o << "\\n"; break;
                case '\r': o << "\\r"; break;
                case '\t': o << "\\t"; break;
                default:
                    if ('\x00' <= c && c <= '\x1f') {
                        o << "\\u" << hex << setw(4) << setfill('0') << (int)c;
                    } else {
                        o << c;
                    }
            }
        }
        return o.str();
    }
    
    /**
     * Writes the snarl tree map to a JSON file.
     */
    void write_json(const string& output_path) {
        ofstream out(output_path);
        if (!out.is_open()) {
            throw runtime_error("Cannot open output file: " + output_path);
        }
        
        out << "{\n";
        bool first = true;
        for (const auto& [key, entry] : snarl_tree_map) {
            if (!first) {
                out << ",\n";
            }
            first = false;
            
            out << "    \"" << json_escape(key) << "\": [\n";
            out << "        [\n";
            
            // Write children IDs
            for (size_t i = 0; i < entry.children_ids.size(); ++i) {
                out << "            \"" << json_escape(entry.children_ids[i]) << "\"";
                if (i < entry.children_ids.size() - 1) {
                    out << ",";
                }
                out << "\n";
            }
            
            out << "        ],\n";
            out << "        " << entry.leaf_count << "\n";
            out << "    ]";
        }
        out << "\n}\n";
        
        out.close();
    }
    
    /**
     * Subdivides a large chain into smaller chunks based on leaf snarl count.
     */
    void subdivide_chain(const net_handle_t& chain, const string& chain_path, const string& root_chain_id) {
        auto [chain_start, chain_end, chain_start_rev, chain_end_rev] = get_boundary_nodes_with_orientation(chain);
        
        // Get ordered children of this chain
        vector<net_handle_t> children;
        index.for_each_child(chain, [&](const net_handle_t& child) {
            children.push_back(child);
        });
        
        if (children.empty()) {
            // No children - shouldn't happen for chains, but handle it
            string chain_id = get_id(chain);
            int leaf_count = 0;
            auto it = snarl_tree_map.find(chain_id);
            if (it != snarl_tree_map.end()) {
                leaf_count = it->second.leaf_count;
            }
            chunk_points.push_back({chain_start, chain_end, leaf_count, chain_path, chain_start_rev, chain_end_rev, root_chain_id});
            return;
        }
        
        // Collect chunks locally first, then merge small last chunk if needed
        vector<ChunkPoint> local_chunks;
        
        nid_t chunk_start_node = chain_start;
        bool chunk_start_is_reverse = chain_start_rev;
        int cumulative_leaf_snarls = 0;
        
        for (size_t i = 0; i < children.size(); ++i) {
            const auto& child = children[i];
            string child_id = get_id(child);
            
            // Get leaf snarl count for this child
            int child_leaf_count = 0;
            auto it = snarl_tree_map.find(child_id);
            if (it != snarl_tree_map.end()) {
                child_leaf_count = it->second.leaf_count;
            }
            
            cumulative_leaf_snarls += child_leaf_count;
            
            // Only create chunk boundaries at chains (nodes), not at snarls
            bool is_child_chain = index.is_chain(child) || index.is_node(child);
            bool reached_target = cumulative_leaf_snarls >= target_leaf_snarls_per_chunk;
            bool is_last_child = (i == children.size() - 1);
            
            // Create chunk boundary only if:
            // 1. We're at a chain/node AND reached target, OR
            // 2. It's the last child (must close the chunk)
            if ((is_child_chain && reached_target) || is_last_child) {
                // Get the end boundary for this chunk
                nid_t chunk_end_node;
                bool chunk_end_is_reverse;
                
                if (index.is_node(child)) {
                    // Child is a node directly
                    handle_t child_handle = index.get_handle(child, &graph);
                    chunk_end_node = graph.get_id(child_handle);
                    chunk_end_is_reverse = graph.get_is_reverse(child_handle);
                } else if (index.is_chain(child)) {
                    // Child is a chain - use its end boundary
                    auto [child_start, child_end, child_start_rev, child_end_rev] = get_boundary_nodes_with_orientation(child);
                    chunk_end_node = child_end;
                    chunk_end_is_reverse = child_end_rev;
                } else if (index.is_snarl(child)) {
                    // Last child is a snarl - use its end boundary
                    auto [child_start, child_end, child_start_rev, child_end_rev] = get_boundary_nodes_with_orientation(child);
                    chunk_end_node = child_end;
                    chunk_end_is_reverse = child_end_rev;
                } else {
                    // Skip unknown types
                    continue;
                }
                
                local_chunks.push_back({chunk_start_node, chunk_end_node, cumulative_leaf_snarls, chain_path, chunk_start_is_reverse, chunk_end_is_reverse, root_chain_id});
                
                // Start new chunk if not the last child
                if (!is_last_child && reached_target) {
                    chunk_start_node = chunk_end_node;
                    chunk_start_is_reverse = chunk_end_is_reverse;
                    cumulative_leaf_snarls = 0;
                }
            }
        }
        
        // Merge last chunk with previous if it's too small (< 30000 leaf snarls)
        const int MIN_LAST_CHUNK_SIZE = 35000;
        if (local_chunks.size() > 1) {
            ChunkPoint& last_chunk = local_chunks.back();
            if (last_chunk.leaf_count < MIN_LAST_CHUNK_SIZE) {
                // Merge with previous chunk
                ChunkPoint& prev_chunk = local_chunks[local_chunks.size() - 2];
                prev_chunk.end_node = last_chunk.end_node;
                prev_chunk.end_is_reverse = last_chunk.end_is_reverse;
                prev_chunk.leaf_count += last_chunk.leaf_count;
                local_chunks.pop_back();  // Remove the merged chunk
            }
        }
        
        // Add all chunks to the global list
        chunk_points.insert(chunk_points.end(), local_chunks.begin(), local_chunks.end());
    }
    
public:
    /**
     * Constructor.
     */
    SnarlTreeBuilder() : target_leaf_snarls_per_chunk(100000) {}
    
    /**
     * Sets the target number of leaf snarls per chunk for chunk point generation.
     */
    void set_target_leaf_snarls_per_chunk(int target) {
        target_leaf_snarls_per_chunk = target;
    }

    /**
     * Loads the graph and snarl distance index from disk.
     */
    void load(const string& graph_path, const string& index_path) {
        auto t0 = chrono::high_resolution_clock::now();
        
        ifstream graph_file(graph_path);
        if (!graph_file.is_open()) {
            throw runtime_error("Cannot open graph file: " + graph_path);
        }
        graph.deserialize(graph_file);
        graph_file.close();
        
        auto t1 = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = t1 - t0;
        cerr << "Loading graph: " << fixed << setprecision(2) 
             << elapsed.count() << " seconds" << endl;
        
        ifstream index_file(index_path);
        if (!index_file.is_open()) {
            throw runtime_error("Cannot open index file: " + index_path);
        }
        index.deserialize(index_file);
        index_file.close();
        
        auto t2 = chrono::high_resolution_clock::now();
        elapsed = t2 - t1;
        cerr << "Loading index: " << fixed << setprecision(2) 
             << elapsed.count() << " seconds" << endl;
    }
    
    /**
     * Loads the snarl tree map from an existing JSON file.
     */
    void load_json(const string& input_path) {
        auto start = chrono::high_resolution_clock::now();
        
        ifstream in(input_path);
        if (!in.is_open()) {
            throw runtime_error("Cannot open input JSON file: " + input_path);
        }
        
        snarl_tree_map.clear();
        
        string line;
        string current_key;
        vector<string> children_ids;
        int leaf_count = 0;
        int state = 0; // 0=looking for key, 1=reading children, 2=reading leaf_count
        
        while (getline(in, line)) {
            // Skip empty lines and braces
            if (line.empty() || line.find_first_not_of(" \t") == string::npos) continue;
            if (line.find("{") != string::npos && line.find("\"") == string::npos) continue;
            if (line.find("}") != string::npos && line.find("\"") == string::npos) continue;
            
            // Parse key
            size_t key_start = line.find("\"");
            if (key_start != string::npos && state == 0) {
                size_t key_end = line.find("\"", key_start + 1);
                current_key = line.substr(key_start + 1, key_end - key_start - 1);
                state = 1;
                children_ids.clear();
                continue;
            }
            
            // Parse children
            if (state == 1) {
                size_t child_start = line.find("\"");
                if (child_start != string::npos) {
                    size_t child_end = line.find("\"", child_start + 1);
                    string child_id = line.substr(child_start + 1, child_end - child_start - 1);
                    children_ids.push_back(child_id);
                } else if (line.find("]") != string::npos) {
                    state = 2; // Next line should be leaf count
                }
                continue;
            }
            
            // Parse leaf count
            if (state == 2) {
                // Find the number in the line
                size_t num_start = line.find_first_of("0123456789");
                if (num_start != string::npos) {
                    size_t num_end = line.find_first_not_of("0123456789", num_start);
                    string num_str = line.substr(num_start, num_end - num_start);
                    leaf_count = stoi(num_str);
                    
                    // Save entry
                    snarl_tree_map[current_key] = {children_ids, leaf_count};
                    state = 0;
                }
            }
        }
        
        in.close();
        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        cerr << "Loaded snarl tree from JSON: " << fixed << setprecision(2) 
             << elapsed.count() << " seconds (" << snarl_tree_map.size() << " entries)" << endl;
    }
    
    /**
     * Builds the snarl tree map using parallel threads.
     */
    void build_parallel(const string& output_path, int num_threads_override = 0) {
        auto workflow_start = chrono::high_resolution_clock::now();
        
        vector<net_handle_t> handles_to_process = collect_handles_at_depth(PARALLELIZATION_DEPTH);
        
        unsigned int num_threads;
        if (num_threads_override > 0) {
            num_threads = num_threads_override;
        } else {
            num_threads = min((unsigned int)thread::hardware_concurrency(), 
                             (unsigned int)handles_to_process.size());
        }
        if (num_threads == 0) num_threads = 1;
        
        cerr << "Running with " << num_threads << " parallel threads." << endl;
        
        if (!handles_to_process.empty()) {
            // Use a thread pool approach to avoid creating too many threads
            vector<thread> threads;
            size_t handles_per_thread = (handles_to_process.size() + num_threads - 1) / num_threads;
            
            for (size_t t = 0; t < num_threads; ++t) {
                threads.emplace_back([this, t, &handles_to_process, handles_per_thread, num_threads]() {
                    size_t start_idx = t * handles_per_thread;
                    size_t end_idx = min(start_idx + handles_per_thread, handles_to_process.size());
                    
                    for (size_t i = start_idx; i < end_idx; ++i) {
                        process_subtree_worker(i, handles_to_process[i]);
                    }
                });
            }
            
            // Wait for all threads to complete
            for (auto& thread : threads) {
                thread.join();
            }
        } else {
            // If no handles to process, traverse from root sequentially
            snarl_tree_map = traverse_decomposition(index.get_root());
        }
        
        // Reconstruct the top levels of the tree
        build_top_tree(index.get_root(), 0);
        
        // Write output
        write_json(output_path);
        
        auto workflow_end = chrono::high_resolution_clock::now();
        chrono::duration<double> total_elapsed = workflow_end - workflow_start;
        
        cerr << "Building snarl tree: " << fixed << setprecision(2) 
             << total_elapsed.count() << " seconds (produced " 
             << snarl_tree_map.size() << " entries)" << endl;
    }
    
    /**
     * Generates chunk points for parallel processing.
     * Each chunk contains approximately target_leaf_snarls_per_chunk leaf snarls.
     */
    void generate_chunk_points() {
        chunk_points.clear();
        
        net_handle_t root = index.get_root();
        
        // Get all children of root (the top-level chains)
        vector<net_handle_t> root_chains;
        index.for_each_child(root, [&](const net_handle_t& child) {
            if (!index.is_node(child)) {
                root_chains.push_back(child);
            }
        });
        
        cerr << "Processing " << root_chains.size() << " root-level chains for chunking..." << endl;
        
        for (const auto& chain : root_chains) {
            string chain_id = get_id(chain);
            
            // Get leaf snarl count for this chain
            auto it = snarl_tree_map.find(chain_id);
            if (it == snarl_tree_map.end()) {
                cerr << "Warning: Chain " << chain_id << " not found in snarl_tree_map" << endl;
                continue;
            }
            
            int chain_leaf_count = it->second.leaf_count;
            
            // Skip chains with no leaf snarls
            if (chain_leaf_count == 0) {
                continue;
            }
            
            // Get the path (CHM13 reference contig) for this chain
            // Use the start boundary node to determine the path
            auto [start_node, end_node, start_is_reverse, end_is_reverse] = get_boundary_nodes_with_orientation(chain);
            string chain_path = get_path_for_node(start_node);
            
            if (chain_leaf_count <= target_leaf_snarls_per_chunk) {
                // Chain is small enough - add it as a single chunk
                chunk_points.push_back({start_node, end_node, chain_leaf_count, chain_path, start_is_reverse, end_is_reverse, chain_id});
            } else {
                // Chain needs to be subdivided
                subdivide_chain(chain, chain_path, chain_id);
            }
        }
        
        cerr << "Generated " << chunk_points.size() << " chunk points (target: ~" 
             << target_leaf_snarls_per_chunk << " leaf snarls per chunk)" << endl;
    }
    
    /**
     * Writes chunk points to a file.
     */
    void write_chunk_points(const string& output_path) {
        ofstream out(output_path);
        if (!out.is_open()) {
            throw runtime_error("Cannot open chunk points output file: " + output_path);
        }
        
        // Write as TSV (tab-separated values)
        out << "start_node\tstart_orientation\tend_node\tend_orientation\tleaf_snarls\tpath_name\troot_chain_id\n";
        for (const auto& chunk : chunk_points) {
            out << chunk.start_node << "\t" 
                << (chunk.start_is_reverse ? "-" : "+") << "\t"
                << chunk.end_node << "\t" 
                << (chunk.end_is_reverse ? "-" : "+") << "\t"
                << chunk.leaf_count << "\t" 
                << chunk.path_name << "\t"
                << chunk.root_chain_id << "\n";
        }
        
        out.close();
        cerr << "Chunk points written to: " << output_path << endl;
    }
};

int main(int argc, char* argv[]) {
    // Simple argument parsing
    string graph_path, index_path, output_path = "snarl_tree_map.json";
    string input_json_path = "";
    string chunk_points_path = "";
    int num_threads = 0;  // 0 means use default (hardware_concurrency)
    int target_chunk_size = 100000;  // Default 100k leaf snarls per chunk
    bool generate_chunks = false;
    
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if ((arg == "-g" || arg == "--graph") && i + 1 < argc) {
            graph_path = argv[++i];
        } else if ((arg == "-i" || arg == "--index") && i + 1 < argc) {
            index_path = argv[++i];
        } else if ((arg == "-j" || arg == "--input-json") && i + 1 < argc) {
            input_json_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output-json") && i + 1 < argc) {
            output_path = argv[++i];
        } else if ((arg == "-n" || arg == "--num-threads") && i + 1 < argc) {
            num_threads = stoi(argv[++i]);
        } else if ((arg == "-c" || arg == "--chunk-points") && i + 1 < argc) {
            chunk_points_path = argv[++i];
            generate_chunks = true;
        } else if ((arg == "-t" || arg == "--target-chunk-size") && i + 1 < argc) {
            target_chunk_size = stoi(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            cout << "Usage: " << argv[0] << " [OPTIONS]\n";
            cout << "\nMode 1: Build snarl tree from graph\n";
            cout << "  Required: -g GRAPH -i INDEX\n";
            cout << "  Optional: -o OUTPUT -n THREADS\n";
            cout << "\nMode 2: Load existing snarl tree and generate chunk points\n";
            cout << "  Required: -j INPUT_JSON -g GRAPH -i INDEX -c CHUNK_FILE\n";
            cout << "  Optional: -t TARGET_SIZE\n";
            cout << "\nMode 3: Build snarl tree and generate chunk points in one step\n";
            cout << "  Required: -g GRAPH -i INDEX -c CHUNK_FILE\n";
            cout << "  Optional: -o OUTPUT -n THREADS -t TARGET_SIZE\n";
            cout << "\nArguments:\n";
            cout << "  -g, --graph              Path to the variation graph (.pg file)\n";
            cout << "  -i, --index              Path to the snarl distance index (.dist file)\n";
            cout << "  -j, --input-json         Load existing snarl tree JSON (skip building)\n";
            cout << "  -o, --output-json        Path to output snarl tree JSON (default: snarl_tree_map.json)\n";
            cout << "  -n, --num-threads        Number of threads to use (default: auto-detect)\n";
            cout << "  -c, --chunk-points       Path to output chunk points TSV file\n";
            cout << "  -t, --target-chunk-size  Target leaf snarls per chunk (default: 100000)\n";
            cout << "  -h, --help               Show this help message\n";
            cout << "\nExamples:\n";
            cout << "  # Build snarl tree only\n";
            cout << "  " << argv[0] << " -g graph.pg -i index.dist -o tree.json\n";
            cout << "\n  # Load existing tree and generate chunks (fast!)\n";
            cout << "  " << argv[0] << " -j tree.json -g graph.pg -i index.dist -c chunks.tsv\n";
            cout << "\n  # Build tree and generate chunks in one step\n";
            cout << "  " << argv[0] << " -g graph.pg -i index.dist -c chunks.tsv\n";
            return 0;
        }
    }
    
    // Validate arguments
    if (!input_json_path.empty()) {
        // Mode 2: Loading existing JSON for chunk generation
        if (graph_path.empty() || index_path.empty()) {
            cerr << "Error: When using -j/--input-json, you must also provide -g/--graph and -i/--index\n";
            cerr << "       (needed for extracting boundary nodes for chunks)\n";
            cerr << "Use -h or --help for usage information.\n";
            return 1;
        }
    } else {
        // Mode 1 or 3: Building snarl tree
        if (graph_path.empty() || index_path.empty()) {
            cerr << "Error: Both -g/--graph and -i/--index are required.\n";
            cerr << "Use -h or --help for usage information.\n";
            return 1;
        }
    }
    
    try {
        SnarlTreeBuilder decomp;
        decomp.set_target_leaf_snarls_per_chunk(target_chunk_size);
        
        // Always load graph and index (needed for chunk boundary extraction)
        decomp.load(graph_path, index_path);
        
        // Either load existing JSON or build new tree
        if (!input_json_path.empty()) {
            decomp.load_json(input_json_path);
        } else {
            decomp.build_parallel(output_path, num_threads);
        }
        
        // Generate chunk points if requested
        if (generate_chunks) {
            decomp.generate_chunk_points();
            decomp.write_chunk_points(chunk_points_path);
        }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}

