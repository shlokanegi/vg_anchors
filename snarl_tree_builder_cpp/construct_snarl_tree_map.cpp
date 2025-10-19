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

class SnarlTreeBuilder {
private:
    PackedGraph graph;
    SnarlDistanceIndex index;
    map<string, TreeMapEntry> snarl_tree_map;
    mutex map_mutex;  // For thread-safe map operations
    
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
    
public:
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
     * Builds the snarl tree map using parallel threads.
     */
    void build_parallel(const string& output_path) {
        auto workflow_start = chrono::high_resolution_clock::now();
        
        vector<net_handle_t> handles_to_process = collect_handles_at_depth(PARALLELIZATION_DEPTH);
        
        unsigned int num_threads = min((unsigned int)thread::hardware_concurrency(), 
                                      (unsigned int)handles_to_process.size());
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
};

int main(int argc, char* argv[]) {
    // Simple argument parsing
    string graph_path, index_path, output_path = "snarl_tree_map.json";
    
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if ((arg == "-g" || arg == "--graph") && i + 1 < argc) {
            graph_path = argv[++i];
        } else if ((arg == "-i" || arg == "--index") && i + 1 < argc) {
            index_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output-json") && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            cout << "Usage: " << argv[0] << " -g GRAPH -i INDEX [-o OUTPUT]\n";
            cout << "  -g, --graph       Path to the variation graph (.pg file)\n";
            cout << "  -i, --index       Path to the snarl distance index (.dist file)\n";
            cout << "  -o, --output-json Path to the output JSON file (default: snarl_tree_map.json)\n";
            return 0;
        }
    }
    
    if (graph_path.empty() || index_path.empty()) {
        cerr << "Error: Both -g/--graph and -i/--index are required.\n";
        cerr << "Use -h or --help for usage information.\n";
        return 1;
    }
    
    try {
        SnarlTreeBuilder decomp;
        decomp.load(graph_path, index_path);
        decomp.build_parallel(output_path);
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}

