//------------------------------------------------------------------------------
// contig_filter_extract.cpp
//
// This program extracts contigs from a GFA assembly file based on minimizer
// hit counts. It identifies outlier contigs (those with high minimizer hit
// counts) by calculating percentiles from the distribution of hits, then
// outputs the selected contigs as FASTA sequences along with metadata and
// visualization tools.
//
// Input: contig_minimizers.tsv (minimizer statistics), GFA assembly file
// Output: FASTA sequences, TSV metadata, R script for visualization
//------------------------------------------------------------------------------

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>
#include <queue>

//------------------------------------------------------------------------------
// Data Structures
//------------------------------------------------------------------------------

/**
 * Metadata for a contig parsed from contig_minimizers.tsv.
 * 
 * This structure stores information about minimizer hits on each contig:
 * - unique_minimizers: Number of distinct minimizers that hit this contig
 * - min_offset/max_offset: Range of positions where minimizers were found
 * 
 * Used to filter contigs based on their minimizer hit frequency.
 */
struct ContigInfo
{
  std::string contig_name;        // Original contig identifier from GFA
  std::size_t unique_minimizers;  // Count of unique minimizers hitting this contig
  std::size_t min_offset;         // Minimum position where a minimizer was found
  std::size_t max_offset;         // Maximum position where a minimizer was found
};

/**
 * Edge representation in the contig adjacency graph parsed from GFA L-lines.
 * 
 * GFA format uses 'L' (link) lines to represent overlaps between contigs:
 *   L <source> <source_orientation> <sink> <sink_orientation> <overlap_length>
 * 
 * This structure stores:
 * - source_contig_name: The source contig in the overlap
 * - sink_contig_name: The destination/target contig in the overlap
 * - overlap_length: Length of the overlap between the two contigs
 * - source_orientation: Orientation of source contig ('+' or '-')
 * - sink_orientation: Orientation of sink contig ('+' or '-')
 * 
 * Used to build a graph representation of contig relationships for distance-based filtering.
 */
struct ContigEdge
{
  std::string source_contig_name;
  std::string sink_contig_name;
  std::size_t overlap_length;
  char source_orientation;
  char sink_orientation;
};
/**
 * A labeled contig with its sequence and metadata for output.
 * 
 * Contains the full sequence along with labeling information:
 * - label: New label format: C<chunk_id>#<contig_name> (e.g., C0#0-17-4-0-P1)
 * - contig_name: Original contig name from the GFA file
 * - sequence: Full nucleotide sequence
 * - original_len/output_len: Sequence length information (currently same since no clipping)
 */
struct LabeledContig
{
  std::string label;        // New label: C<chunk_id>#<contig_name>
  std::string contig_name;  // Original contig name from GFA
  std::string sequence;     // Full sequence (or clipped if requested in future)
  std::size_t original_len; // Original length before any clipping
  std::size_t output_len;   // Length of output sequence
};

/**
 * A path sequence representing a concatenated and clipped candidate path.
 * 
 * Contains the final sequence for a candidate path after:
 * - Concatenating contig sequences while accounting for overlaps
 * - Clipping the first and last contigs based on minimizer offsets
 * - Keeping middle contigs entirely
 */
struct PathSequence
{
  std::string label;           // Path label: C<chunk_id>#Path<path_id>:<contigs>
  std::vector<std::string> path_contigs;  // Ordered list of contig names in the path
  std::string sequence;        // Concatenated and clipped sequence
  std::size_t original_length; // Length before clipping (concatenated length)
  std::size_t output_length;   // Length after clipping
  std::vector<std::size_t> original_contig_lengths;  // Original length of each contig in path
  std::vector<std::size_t> clipped_contig_lengths;   // Clipped length of each contig in path (first/last clipped, middle full)
  std::size_t min_offset_in_start_node;  // min_offset of the first contig in path
  std::size_t max_offset_in_sink_node;   // max_offset of the last contig in path
};


//------------------------------------------------------------------------------
// Parsing Functions
//------------------------------------------------------------------------------

/**
 * Parse a comma-separated list of integers (offset positions).
 * 
 * Used to extract minimizer offset positions from the TSV file format.
 * Example input: "10,25,47,92" -> [10, 25, 47, 92]
 * 
 * @param offsets_str Comma-separated string of integers
 * @return Vector of parsed offset values
 */
std::vector<std::size_t> parse_offsets(const std::string& offsets_str)
{
  std::vector<std::size_t> result;
  if(offsets_str.empty()) { return result; }
  
  // Tokenize by comma delimiter
  std::stringstream ss(offsets_str);
  std::string token;
  while(std::getline(ss, token, ','))
  {
    if(!token.empty())
    {
      // Convert string to unsigned long long, then cast to size_t
      result.push_back(static_cast<std::size_t>(std::stoull(token)));
    }
  }
  return result;
}

/**
 * Load contig_minimizers.tsv file into a map keyed by contig name.
 * 
 * Expected TSV format (tab-separated):
 *   contig_name <tab> unique_minimizer_count <tab> minimizer_indices <tab> contig_offsets
 * 
 * The function parses the offset column (comma-separated integers) and computes
 * the min/max offset range for each contig. This helps identify where on the
 * contig the minimizers are located.
 * 
 * @param path Path to the contig_minimizers.tsv file
 * @return Map from contig name to ContigInfo metadata
 * @throws std::runtime_error if file cannot be opened
 */
std::unordered_map<std::string, ContigInfo>
load_contig_filter(const std::string& path)
{
  std::unordered_map<std::string, ContigInfo> result;
  std::ifstream input(path);
  if(!input)
  {
    throw std::runtime_error("Could not open " + path);
  }
  
  // Skip header line
  std::string header;
  if(!std::getline(input, header)) { return result; }

  // Parse each data line
  std::string line;
  while(std::getline(input, line))
  {
    if(line.empty()) { continue; }
    
    // Tokenize by tab delimiter
    std::stringstream ss(line);
    std::string contig, count_str, minimizer_indices, contig_offsets;
    
    // Extract the four tab-separated fields
    if(!std::getline(ss, contig, '\t')) { continue; }
    if(!std::getline(ss, count_str, '\t')) { continue; }
    if(!std::getline(ss, minimizer_indices, '\t')) { continue; }
    if(!std::getline(ss, contig_offsets, '\t')) { continue; }
    
    // Parse minimizer count and offset positions
    std::size_t count = std::stoull(count_str);
    std::vector<std::size_t> offsets = parse_offsets(contig_offsets);
    
    // Calculate offset range (min and max positions where minimizers were found)
    std::size_t min_offset = SIZE_MAX;
    std::size_t max_offset = 0;
    for(std::size_t off : offsets)
    {
      if(off < min_offset) { min_offset = off; }
      if(off > max_offset) { max_offset = off; }
    }
    // Handle case where no offsets were found
    if(min_offset == SIZE_MAX) { min_offset = 0; }
    
    // Store contig metadata
    result[contig] = { contig, count, min_offset, max_offset };
  }
  return result;
}

/**
 * Load GFA (Graphical Fragment Assembly) sequences from S-lines.
 * 
 * GFA format uses 'S' lines to define segments (contigs):
 *   S <segment_name> <sequence> <optional_fields>
 * 
 * This function extracts only the segment names and their nucleotide sequences,
 * ignoring other GFA line types (L, H, P, etc.) and optional fields.
 * 
 * @param gfa_path Path to the GFA assembly file
 * @return Map from segment/contig name to nucleotide sequence
 * @throws std::runtime_error if file cannot be opened
 */
std::unordered_map<std::string, std::string>
load_gfa_sequences(const std::string& gfa_path)
{
  std::unordered_map<std::string, std::string> sequences;
  std::ifstream input(gfa_path);
  if(!input)
  {
    throw std::runtime_error("Could not open " + gfa_path);
  }

  std::string line;
  while(std::getline(input, line))
  {
    // Skip non-S lines (S = segment/sequence line in GFA format)
    if(line.empty() || line[0] != 'S') { continue; }
    
    // Parse S-line: S <name> <sequence> [optional fields]
    std::stringstream ss(line);
    std::string type, name, sequence;
    if(!std::getline(ss, type, '\t')) { continue; }  // 'S'
    if(!std::getline(ss, name, '\t')) { continue; }   // segment name
    if(!std::getline(ss, sequence, '\t')) { continue; } // sequence
    
    sequences[name] = sequence;
  }
  return sequences;
}

//------------------------------------------------------------------------------
// Distribution Utilities
//------------------------------------------------------------------------------

/**
 * Calculate the percentile value from a sorted vector using linear interpolation.
 * 
 * Uses the standard linear interpolation method (Type 7) to calculate percentiles.
 * For example, the 90th percentile means 90% of values are below this threshold.
 * 
 * Algorithm:
 *   1. Calculate index position: index = (percentile/100) * (n-1)
 *   2. Find bounding values at floor(index) and ceil(index)
 *   3. Interpolate between them based on fractional part
 * 
 * @param sorted_values Sorted vector of values (must be in ascending order, non-empty)
 * @param percentile Percentile to calculate (0.0-100.0)
 * @return Value at the specified percentile
 * @throws std::runtime_error if input vector is empty
 */
std::size_t calculate_percentile(const std::vector<std::size_t>& sorted_values, double percentile)
{
  if(sorted_values.empty())
  {
    throw std::runtime_error("Cannot calculate percentile of empty distribution");
  }
  
  // Handle edge cases: 0th percentile = minimum, 100th percentile = maximum
  if(percentile <= 0.0) { return sorted_values.front(); }
  if(percentile >= 100.0) { return sorted_values.back(); }
  
  // Linear interpolation method for percentile calculation
  // Formula: index = (percentile / 100) * (n - 1)
  // This gives us the position in the sorted array where the percentile value lies
  double index = (percentile / 100.0) * (sorted_values.size() - 1);
  std::size_t lower_idx = static_cast<std::size_t>(index);
  std::size_t upper_idx = lower_idx + 1;
  
  // Guard against index out of bounds
  if(upper_idx >= sorted_values.size())
  {
    return sorted_values.back();
  }
  
  // Interpolate between lower and upper values
  // weight = fractional part of index (0.0 to 1.0)
  double weight = index - lower_idx;
  return static_cast<std::size_t>(
    sorted_values[lower_idx] * (1.0 - weight) + 
    sorted_values[upper_idx] * weight
  );
}

/**
 * Calculate IQR-based outlier threshold using the interquartile range method.
 * 
 * IQR (Interquartile Range) method identifies outliers by:
 *   1. Calculate Q1 (25th percentile) and Q3 (75th percentile)
 *   2. Calculate IQR = Q3 - Q1
 *   3. Upper fence = Q3 + (multiplier * IQR)
 * 
 * Contigs with values above the upper fence are considered outliers.
 * Common multipliers:
 *   - 1.5: Mild outliers (standard)
 *   - 3.0: Extreme outliers (more stringent)
 * 
 * @param sorted_values Sorted vector of values (must be in ascending order, non-empty)
 * @param multiplier Multiplier for IQR (e.g., 1.5 for mild, 3.0 for extreme outliers)
 * @return Threshold value above which contigs are considered outliers
 * @throws std::runtime_error if input vector is empty
 */
std::size_t calculate_iqr_threshold(const std::vector<std::size_t>& sorted_values, double multiplier)
{
  if(sorted_values.empty())
  {
    throw std::runtime_error("Cannot calculate IQR threshold of empty distribution");
  }
  
  // Calculate Q1 (25th percentile) and Q3 (75th percentile)
  double q1 = calculate_percentile(sorted_values, 25.0);
  double q3 = calculate_percentile(sorted_values, 75.0);
  
  // Calculate IQR (Interquartile Range)
  double iqr = static_cast<double>(q3) - static_cast<double>(q1);
  
  // Calculate upper fence: Q3 + (multiplier * IQR)
  double threshold = static_cast<double>(q3) + (multiplier * iqr);
  
  // Return as size_t (round up to ensure we don't miss outliers due to truncation)
  return static_cast<std::size_t>(std::ceil(threshold));
}

//------------------------------------------------------------------------------
// Sequence Utilities
//------------------------------------------------------------------------------

/**
 * Compute the reverse complement of a DNA sequence.
 * 
 * Transforms A->T, C->G, G->C, T->A and reverses the sequence.
 * Non-standard nucleotides (N, etc.) are converted to 'N'.
 * Handles both uppercase and lowercase input.
 * 
 * @param seq Input DNA sequence (forward strand)
 * @return Reverse complement sequence
 */
std::string reverse_complement(const std::string& seq)
{
  std::string result;
  result.reserve(seq.size());
  for(auto iter = seq.rbegin(); iter != seq.rend(); ++iter)
  {
    char c = *iter;
    switch(c)
    {
    case 'A': result.push_back('T'); break;
    case 'C': result.push_back('G'); break;
    case 'G': result.push_back('C'); break;
    case 'T': result.push_back('A'); break;
    case 'a': result.push_back('t'); break;
    case 'c': result.push_back('g'); break;
    case 'g': result.push_back('c'); break;
    case 't': result.push_back('a'); break;
    default: result.push_back('N'); break;
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Build Candidate Contig Records
//------------------------------------------------------------------------------

/**
 * Build labeled contig records for contigs that pass the minimizer threshold filter.
 * 
 * This function:
 *   1. Filters contigs based on minimizer hit count (>= threshold)
 *   2. Retrieves their sequences from the GFA data
 *   3. Creates labeled records with format: C<chunk_id>#<contig_name>:<min_offset>-<max_offset>
 *      (or C<chunk_id>#<contig_name> if keeping full sequence)
 *   4. Chunks sequences to include the region between min_offset and max_offset
 *      with 50bp flanking on both sides, unless the clipped length is >= 
 *      keep_full_threshold of the original length, in which case the full contig is kept
 * 
 * Contigs are selected as outliers if their unique_minimizers count meets or
 * exceeds the threshold (typically the 90th percentile or higher). By default,
 * sequences are extracted from (min_offset - 50bp) to (max_offset + 50bp) to
 * include the region where minimizers were found with flanking sequence. However,
 * if the clipped length would be >= keep_full_threshold (default: 90%) of the 
 * original contig length, the full contig is retained to avoid unnecessary fragmentation.
 * 
 * @param contigs Map of all contig metadata (minimizer info)
 * @param sequences Map of contig name to sequence from GFA file
 * @param chunk_id Chunk identifier for labeling (e.g., "0", "1")
 * @param threshold Minimum minimizer hit count required (percentile-based)
 * @param keep_full_threshold Fraction threshold (0.0-1.0) for keeping full contig (default: 0.90)
 *                            If clipped length >= keep_full_threshold * original_length, keep full
 * @return Vector of labeled contig records, sorted by contig name
 */
/**
 * Filter candidate contigs based on minimizer hit count threshold.
 * 
 * This function performs the first filtering step: selecting contigs whose
 * unique_minimizers count meets or exceeds the specified threshold. The threshold
 * is typically calculated from a percentile (e.g., 90th percentile) or IQR-based
 * outlier detection method.
 * 
 * This is the initial filtering step before distance-based filtering. Only contigs
 * that pass this threshold filter are considered for further analysis in the
 * contig adjacency graph.
 * 
 * @param contigs Map of all contig metadata (minimizer info) from TSV file
 * @param threshold Minimum minimizer hit count required to be a candidate
 *                 (calculated from percentile or IQR method)
 * @param keep_full_threshold Unused in this function (kept for API consistency)
 * @return Vector of ContigInfo records that passed the threshold filter
 */
std::vector<ContigInfo>
filter_candidate_contig_records_by_minimizer_threshold(const std::unordered_map<std::string, ContigInfo>& contigs,
                     std::size_t threshold,
                     double keep_full_threshold = 0.90)
{
  std::vector<ContigInfo> candidate_contigs;
  
  // Iterate through all contigs and filter by threshold
  for(const auto& item : contigs)
  {
    const ContigInfo& info = item.second;
    
    // Filter: skip contigs with minimizer counts below the threshold
    // Threshold is typically calculated from percentile (e.g., 90th percentile)
    if(info.unique_minimizers < threshold) { continue; }
    candidate_contigs.push_back(info);
  }
  // If no candidate contigs found, then threshold is too high. So, skip this approach, and pick all contigs as candidate contigs
  if (candidate_contigs.size() == 0) {
    std::cout << "No candidate contigs found after filtering by minimizer threshold. Picking all contigs as candidate contigs" << std::endl;
    std::vector<ContigInfo> contig_info_list;
    for (const auto& item : contigs) {
      const ContigInfo& info = item.second;
      contig_info_list.push_back(info);
    }
    std::sort(contig_info_list.begin(), contig_info_list.end(), [](const ContigInfo& a, const ContigInfo& b) {
      return a.unique_minimizers > b.unique_minimizers;
    });
    for (int i = 0; i < contig_info_list.size(); i++) {
      candidate_contigs.push_back(contig_info_list[i]);
    }
  }
  // Debug output: print filtered candidate contigs
  std::cout << "================================================================" << std::endl;
  std::cout << "FILTERING CANDIDATES" << std::endl;
  std::cout << "================================================================" << std::endl;

  std::cout << "  filtered candidate_contigs by minimizer threshold. Num remaining: " << candidate_contigs.size() << std::endl;
  std::cout << "  candidate_contigs: " << std::endl;
  for (const auto& info : candidate_contigs) {
    std::cout << "  " << info.contig_name << " : " << info.unique_minimizers << std::endl;
  }
  std::cout << std::endl;
  return candidate_contigs;
}

/**
 * Build an adjacency graph representation of contigs from GFA L-lines.
 * 
 * Parses the GFA file and extracts all 'L' (link) lines which represent overlaps
 * between contigs. Builds a directed graph where:
 * - Nodes: contig names
 * - Edges: overlaps between contigs (from L-lines)
 * 
 * GFA L-line format:
 *   L <source> <source_orientation> <sink> <sink_orientation> <overlap_length>
 * Example: L contig1 + contig2 + 1000
 * 
 * The graph is represented as an adjacency list: each source contig maps to a
 * vector of ContigEdge structures representing all outgoing edges from that contig.
 * 
 * This graph is used by select_connected_candidate_contigs to perform distance-based
 * filtering to find connected components of candidate contigs within max_distance hops.
 * 
 * @param contig_gfa_path Path to the GFA assembly file containing L-lines
 * @return Adjacency list: map from source contig name to vector of outgoing edges
 */
std::unordered_map<std::string, std::vector<ContigEdge> >
build_contig_graph(const std::string& contig_gfa_path)
{
  std::unordered_map<std::string, std::vector<ContigEdge> > contig_graph;
  std::ifstream in(contig_gfa_path);
  std::string line;
  
  // Parse each line in the GFA file
  while(std::getline(in, line))
  {
    std::string link_line_symbol;
    if (line[0] == 'S') {
      // extract the contig name from the S-line. And make an entry in the contig_graph with the contig name as the key, and an empty vector of ContigEdge as the value.
      std::string contig_name;
      std::string sequence;
      std::string length_tag;
      std::stringstream ss(line);
      ss >> link_line_symbol >> contig_name >> sequence >> length_tag;
      contig_graph[contig_name] = std::vector<ContigEdge>();
    }
    else if(line[0] == 'L')
    {
      std::string source_contig_name, sink_contig_name;
      char source_orientation, sink_orientation;
      std::size_t overlap_length;
      
      // Parse L-line: L <source> <source_orientation> <sink> <sink_orientation> <overlap_length>
      std::stringstream ss(line);
      ss >> link_line_symbol >> source_contig_name >> source_orientation >> sink_contig_name >> sink_orientation >> overlap_length;
      
      // Initialize adjacency list entry if this is the first edge from source_contig_name
      if (contig_graph.find(source_contig_name) == contig_graph.end()) {
        contig_graph[source_contig_name] = std::vector<ContigEdge>();
      }
      
      // Add edge from source to sink
      contig_graph[source_contig_name].push_back(ContigEdge{source_contig_name, sink_contig_name, overlap_length, source_orientation, sink_orientation});
    }
  }
  return contig_graph;
}

std::vector<std::string> split_string_by_hyphen(const std::string& s) {
  std::vector<std::string> tokens;
  std::string delimiter = "-";
  size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  std::string token;

  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
      token = s.substr(pos_start, pos_end - pos_start);
      pos_start = pos_end + delim_len;
      tokens.push_back(token);
  }

  tokens.push_back(s.substr(pos_start)); // Add the last part
  return tokens;
}

bool check_siblings(const ContigInfo& contig1, const ContigInfo& contig2) {
  std::vector<std::string> contig1_parts = split_string_by_hyphen(contig1.contig_name);
  int n1 = contig1_parts.size();
  std::vector<std::string> contig2_parts = split_string_by_hyphen(contig2.contig_name);
  int n2 = contig2_parts.size();
  
  if (n1 != n2) {
    return false;
  }
  // // let's not enforce that the first parts match (and same for last parts). Because there could be complex examples where this fails.
  // if ((contig1_parts[0] != contig2_parts[0]) || (contig1_parts[n1-1] != contig2_parts[n2-1])) {
  //   return false;
  // }
  int num_different_parts = 0;
  for (int i = 1; i < n1-1; i++) {
    if (contig1_parts[i] != contig2_parts[i]) {
      num_different_parts++;
    }
  }
  return (num_different_parts == 1);
}

/**
 * Convert directed contig graph to bidirected representation.
 * 
 * Takes a directed adjacency graph and creates an undirected (bidirected) version
 * where each edge is represented in both directions. This allows traversal in
 * both directions for finding connected components.
 * 
 * @param directed_contig_graph Directed adjacency list from GFA L-lines
 * @return Bidirected adjacency list with edges in both directions
 */
std::unordered_map<std::string, std::vector<std::string> > 
build_bidirected_contig_graph_from_directed(const std::unordered_map<std::string, std::vector<ContigEdge> >& directed_contig_graph)
{
  std::unordered_map<std::string, std::vector<std::string> > bidirected_contig_graph;
  for (const auto& [contig_name, edges] : directed_contig_graph) {
    if (bidirected_contig_graph.find(contig_name) == bidirected_contig_graph.end()) {
      bidirected_contig_graph[contig_name] = std::vector<std::string>();
    }
    for (const auto& edge : edges) {
      if (bidirected_contig_graph.find(edge.sink_contig_name) == bidirected_contig_graph.end()) {
        bidirected_contig_graph[edge.sink_contig_name] = std::vector<std::string>();
      }
      if (bidirected_contig_graph.find(edge.source_contig_name) == bidirected_contig_graph.end()) {
        bidirected_contig_graph[edge.source_contig_name] = std::vector<std::string>();
      }
      bidirected_contig_graph[edge.sink_contig_name].push_back(edge.source_contig_name);
      bidirected_contig_graph[edge.source_contig_name].push_back(edge.sink_contig_name);
    }
  }
  return bidirected_contig_graph;
}

/**
 * Explore the connected component of a given contig using DFS with distance constraint.
 * 
 * Performs depth-first search to explore a connected component, collecting all
 * candidate contigs reachable within max_distance hops.
 * 
 * @param bidirected_contig_graph Bidirected adjacency list of the contig graph
 * @param curr_contig_name Current contig being explored
 * @param component_contig_names Output vector to collect candidate contigs in this component
 * @param visited_contig_names Map tracking visited contigs to avoid cycles
 * @param candidate_contig_names Set of candidate contig names for quick lookup
 * @param current_distance Current distance from the last candidate contig
 * @param max_distance Maximum distance to explore from candidate contigs
 */
void 
dfs_explore_component_with_given_max_distance(const std::unordered_map<std::string, std::vector<std::string> >& bidirected_contig_graph,
                                                   const std::string& curr_contig_name,
                                                   std::vector<std::string>& component_contig_names,
                                                   std::unordered_map<std::string, bool >& visited_contig_names, 
                                                   std::set<std::string>& candidate_contig_names,
                                                   std::size_t current_distance,
                                                   std::size_t max_distance = 2)
{
  visited_contig_names[curr_contig_name] = true;
  if (candidate_contig_names.find(curr_contig_name) != candidate_contig_names.end()) {
    // add this contig to the component
    std::cout << "Adding contig " << curr_contig_name << " to component" << std::endl;
    component_contig_names.push_back(curr_contig_name);
    // reset the distance counter to 0, to allow the search to continue from this candidate
    current_distance = 0;
  }
  if (current_distance >= max_distance) {
    return;
  }
  for (const auto& neighbor : bidirected_contig_graph.at(curr_contig_name)) {
    if (((visited_contig_names.find(neighbor) != visited_contig_names.end()) && (!visited_contig_names[neighbor])) || (visited_contig_names.find(neighbor) == visited_contig_names.end()))
    {
      dfs_explore_component_with_given_max_distance(bidirected_contig_graph, neighbor, component_contig_names, visited_contig_names, candidate_contig_names, current_distance + 1, max_distance);
    }
  }
}

void topological_sort_helper_dfs_util(const std::unordered_map<std::string, std::vector<ContigEdge> >& contig_graph, const std::string& contig_name, std::vector<std::string>& topological_order, std::unordered_set<std::string>& visited_contig_names, const std::set<std::string>& candidate_contig_names) {
  visited_contig_names.insert(contig_name);
  for (const auto& neighbor : contig_graph.at(contig_name)) {
    if ((visited_contig_names.find(neighbor.sink_contig_name) == visited_contig_names.end()) && (neighbor.sink_contig_name != contig_name)) {
      topological_sort_helper_dfs_util(contig_graph, neighbor.sink_contig_name, topological_order, visited_contig_names, candidate_contig_names);
    }
  }
  auto it = std::find(candidate_contig_names.begin(), candidate_contig_names.end(), contig_name);
  if (it != candidate_contig_names.end()) {
    topological_order.push_back(contig_name);
  }
  return;
}


std::vector<std::string> topological_sort_helper(const std::unordered_map<std::string, std::vector<ContigEdge> >& contig_graph, const std::set<std::string>& candidate_contig_names) {
  std::vector<std::string> topological_order;
  std::unordered_set<std::string> visited_contig_names;
  for (const auto& contig_name : candidate_contig_names) {
    if (visited_contig_names.find(contig_name) == visited_contig_names.end()) {
      topological_sort_helper_dfs_util(contig_graph, contig_name, topological_order, visited_contig_names, candidate_contig_names);
    }
  }
  return topological_order;
}

/**
 * Select connected component of candidate contigs based on distance constraint.
 * 
 * This function finds connected components in the bidirected contig graph where
 * candidate contigs are considered "connected" if they are within max_distance hops
 * of each other (with distance reset to 0 when encountering a candidate contig).
 * It then selects the component containing the candidate contig with the highest
 * unique_minimizers count.
 * 
 * Algorithm:
 *   1. Convert directed graph to bidirected graph to simplify traversal, as connectedness is bidirectional in our definition.
 *   2. For each unvisited candidate contig, perform DFS to explore its connected component
 *   3. During DFS, collect all candidate contigs reachable within max_distance hops
 *      (distance resets to 0 when a candidate contig is encountered)
 *   4. Identify the component containing the candidate contig with highest unique_minimizers
 *   5. Return all candidate contigs from that component
 * 
 * Fallback: If no components are found (no candidates are reachable via graph traversal),
 * falls back to sibling-based grouping using contig name patterns.
 * 
 * Special case: If there's only 1 candidate, it's automatically kept.
 * 
 * @param contig_graph Directed adjacency list representation of contig overlaps (from GFA L-lines)
 * @param candidate_contigs Vector of candidate contigs (already filtered by minimizer threshold)
 * @param max_distance Maximum number of hops allowed between candidate contigs in a component (default: 2)
 * @param contig_type Type of contig ("target" or "query")
 * @param component_minimizer_totals Output parameter: vector of total minimizers per component (tuple format for TSV)
 * @param component_contig_totals Output parameter: vector of number of contigs per component (tuple format for TSV)
 * @return Vector of candidate contigs from the selected connected component
 */
std::vector<ContigInfo>
select_connected_candidate_contigs(const std::unordered_map<std::string, std::vector<ContigEdge> >& contig_graph,
                                 const std::vector<ContigInfo>& candidate_contigs,
                                 std::size_t max_distance,
                                 const std::string& contig_type,
                                 std::vector<std::size_t>& component_minimizer_totals,
                                 std::vector<std::size_t>& component_contig_totals)
{
  std::vector<ContigInfo> output;
  // Build a set of candidate contig names for O(1) lookup
  std::set<std::string> candidate_contig_names;
  for(const auto& info : candidate_contigs) {
    candidate_contig_names.insert(info.contig_name);
  }
  std::unordered_map<std::string, ContigInfo> candidate_contig_info_map;
  for (const auto& info : candidate_contigs) {
    candidate_contig_info_map[info.contig_name] = info;
  }
  
  // Track which candidate contigs are reachable from other candidates
  std::set<std::string> reachable_candidate_contig_names;
  
  // Special case: if only one candidate, keep it automatically
  if (candidate_contigs.size() == 1) {
    reachable_candidate_contig_names.insert(candidate_contigs[0].contig_name);
    // Add total minimizers and contig count for this single-contig component
    component_minimizer_totals.push_back(candidate_contigs[0].unique_minimizers);
    component_contig_totals.push_back(1);
  }
  else {
    std::unordered_map<std::string, bool > visited_contig_names;
    std::unordered_map<std::string, std::vector<std::string> > bidirected_contig_graph = build_bidirected_contig_graph_from_directed(contig_graph);
  
    std::unordered_map<int, std::vector<std::string> > components_of_candidate_contigs;
    int component_idx = 0;
    for (const auto& curr_contig : candidate_contigs) {
      auto curr_contig_name = curr_contig.contig_name;
      if ((visited_contig_names.find(curr_contig_name) != visited_contig_names.end()) && (visited_contig_names[curr_contig_name] == true)) {
        continue;
      }
      std::cout << "Exploring component from contig " << curr_contig_name << std::endl;
      std::vector<std::string> component_contig_names;
      dfs_explore_component_with_given_max_distance(bidirected_contig_graph, curr_contig_name, component_contig_names, visited_contig_names, candidate_contig_names, 0, max_distance);
      std::cout << "Component " << component_idx << " has " << component_contig_names.size() << " contigs" << std::endl;
      components_of_candidate_contigs[component_idx] = component_contig_names;
      
      // Calculate total minimizers for this component
      std::size_t max_minimizers = 0;
      for (const auto& contig_name : component_contig_names) {
        if (candidate_contig_info_map.find(contig_name) != candidate_contig_info_map.end()) {
          max_minimizers = std::max(max_minimizers, candidate_contig_info_map[contig_name].unique_minimizers);
        }
      }
      component_minimizer_totals.push_back(max_minimizers);
      component_contig_totals.push_back(component_contig_names.size());
      
      component_idx++;
    }
  
    // now we add all the contigs of the component that has the highest unique minimizers to the reachable_candidate_contig_names

    /*
    TODO: Change this approach. We shouldn't just consider the highest unique minimizer component. Rather, we should check all components that have close to the highest minimizer count. 
    Then amongst them, we select the component that is closest to the end of the contig graph. 
    The assumption here is that for a component to represent the overlap end, it should be towards the end of the contig graph.
    
    UPDATE: Our simple approach that implements the above idea is the following:- We select the top-2 connected components by highest unique minimizer counts (this top-2 selection is heuristic). 
    Then, if the 2nd component has max unique minimizer count less than half that of the 1st component, we select the 1st component for output (this /2 check is heuristic). 
    (
      This above "divide by 2" thresholding check was added to deal with breakpoints (softclip) in the contig graph, which could cause the wrong component to appear closer to the end of the contig graph.
      We encountered this in chunk 69_70.
    ) 
    Else, among the two, we select the component that is closest to the end of the contig graph.
    */

    std::vector<std::string> topological_order = topological_sort_helper(contig_graph, candidate_contig_names);
    if (contig_type == "query") {
      std::reverse(topological_order.begin(), topological_order.end());
    }
    std::cout<<"Topological order: ";
    for (const auto& contig_name : topological_order) {
      std::cout << contig_name << ";  ";
    }
    std::cout << std::endl;
    
    // we choose top two components by highest unique minimizer count, and then amongst them, we choose the occurs earlier in the topological order.
    std::vector<std::pair<int, int> > components_with_highest_unique_minimizers;  // {unique_minimizers, component_idx}
    for (const auto& [component_idx, component_contig_names] : components_of_candidate_contigs) {
      int highest_unique_minimizers_in_current_component = 0;
      for (const auto& contig_name : component_contig_names) {
        if (candidate_contig_info_map[contig_name].unique_minimizers > highest_unique_minimizers_in_current_component) {  
          highest_unique_minimizers_in_current_component = candidate_contig_info_map[contig_name].unique_minimizers;
        }
      }
      components_with_highest_unique_minimizers.push_back({highest_unique_minimizers_in_current_component, component_idx});
    }
    std::sort(components_with_highest_unique_minimizers.begin(), components_with_highest_unique_minimizers.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
      return a.first > b.first;
    });

    if ((components_with_highest_unique_minimizers.size() < 2) || ((components_with_highest_unique_minimizers[0].first / 2) > components_with_highest_unique_minimizers[1].first)) {
      // pick the component with the highest unique minimizer count
      for (const auto& contig_name : components_of_candidate_contigs[components_with_highest_unique_minimizers[0].second]) {
        reachable_candidate_contig_names.insert(contig_name);
      }
    }
    else {
      // we select the component that is closest to the end of the contig graph.
      std::vector<int>earlier_pos_of_component_in_topological_order;
      for (int idx = 0; idx < components_with_highest_unique_minimizers.size(); idx++) {
        std::vector<std::string> component_contig_names = components_of_candidate_contigs[components_with_highest_unique_minimizers[idx].second];
        int current_component_earliest_pos = INT_MAX;
        for (const auto& contig_name : component_contig_names) {
          auto it = std::find(topological_order.begin(), topological_order.end(), contig_name);
          if (it != topological_order.end()) {
            if (it - topological_order.begin() < current_component_earliest_pos) {
              current_component_earliest_pos = it - topological_order.begin();
            }
          }
        }
        earlier_pos_of_component_in_topological_order.push_back(current_component_earliest_pos);
      }
      std::cout << "Earliest positions of components in topological order: ";
      for (const auto& pos : earlier_pos_of_component_in_topological_order) {
        std::cout << pos << ";  ";
      }
      int selected_component_idx_in_components_with_highest_unique_minimizers_list = 0;
      if (earlier_pos_of_component_in_topological_order.size() >= 2) {
        selected_component_idx_in_components_with_highest_unique_minimizers_list = (earlier_pos_of_component_in_topological_order[0] < earlier_pos_of_component_in_topological_order[1]) ? 0 : 1;
      }
      for (const auto& contig_name : components_of_candidate_contigs[components_with_highest_unique_minimizers[selected_component_idx_in_components_with_highest_unique_minimizers_list].second]) {
        reachable_candidate_contig_names.insert(contig_name);
      }
    }
  }

  // If no components were found, this could mean that the candidate contigs were in sibling groups, and each sibling group doesn't share any links.
  // So, we select the sibling group that has the contig with the highest count of unique minimizers.
  if (reachable_candidate_contig_names.size() <= 1) {
    std::unordered_map<std::string, int> contig_name_to_group_id;
    int incrementing_group_id = 0;
    std::string contig_with_highest_unique_minimizers = candidate_contigs[0].contig_name;
    int highest_unique_minimizers = candidate_contigs[0].unique_minimizers;
    for (const auto& contig : candidate_contigs) {
      if (contig.unique_minimizers > highest_unique_minimizers) {
        highest_unique_minimizers = contig.unique_minimizers;
        contig_with_highest_unique_minimizers = contig.contig_name;
      }
      for (int idx = 0; idx < candidate_contigs.size(); idx++) {
        const ContigInfo& other_contig = candidate_contigs[idx];
        if (other_contig.contig_name == contig.contig_name) {
          break;
        }
        // TODO: Improve sibling checking method. Currently, it only checks contig names for similarity.
        if (check_siblings(contig, other_contig)) {
          contig_name_to_group_id[contig.contig_name] = contig_name_to_group_id[other_contig.contig_name];
          break;
        }
      }
      if (contig_name_to_group_id.find(contig.contig_name) == contig_name_to_group_id.end()) {
        contig_name_to_group_id[contig.contig_name] = ++incrementing_group_id;
      }
    }

    // Calculate total minimizers and contig count for the selected sibling group
    std::size_t sibling_group_total_minimizers = 0;
    std::size_t sibling_group_contig_count = 0;
    for (const auto& contig_name : candidate_contig_names) {
      if (contig_name_to_group_id[contig_name] == contig_name_to_group_id[contig_with_highest_unique_minimizers]) {
        reachable_candidate_contig_names.insert(contig_name);
        sibling_group_contig_count++;
        if (candidate_contig_info_map.find(contig_name) != candidate_contig_info_map.end()) {
          sibling_group_total_minimizers += candidate_contig_info_map[contig_name].unique_minimizers;
        }
      }
    }
    // If we used sibling fallback and component_minimizer_totals is empty, add the sibling group stats
    if (component_minimizer_totals.empty()) {
      component_minimizer_totals.push_back(sibling_group_total_minimizers);
      component_contig_totals.push_back(sibling_group_contig_count);
    }
  }
  // Filter output to only include reachable candidates
  for (const auto& info : candidate_contigs) {
    if (reachable_candidate_contig_names.find(info.contig_name) != reachable_candidate_contig_names.end()) {
      output.push_back(info);
    }
  }  
  // Debug output
  std::cout << "  Kept connected component of candidate contigs with the highest unique minimizers. Num remaining: " << output.size() << std::endl;

  return output;
}

/**
 * Build a subgraph containing only candidate contigs and edges between them.
 * 
 * This function extracts a subgraph from the full contig adjacency graph that
 * only includes:
 *   - Nodes: Only candidate contigs (from filtered_candidate_contigs)
 *   - Edges: Only edges where both source AND sink contigs are candidates
 * 
 * Edges are removed if:
 *   - The source contig is not a candidate (even if sink is a candidate)
 *   - The sink contig is not a candidate (even if source is a candidate)
 * 
 * This creates a "candidate-only" subgraph that will be used for path finding
 * in the next step. The subgraph isolates the candidate contigs and their
 * interconnections, removing connections to non-candidate contigs.
 * 
 * @param full_contig_graph The complete adjacency graph from the GFA file
 * @param candidate_contigs Vector of candidate contigs (after all filtering steps)
 * @return Adjacency list subgraph containing only candidate contigs and edges between them
 */
std::unordered_map<std::string, std::vector<ContigEdge> >
build_candidate_subgraph(const std::unordered_map<std::string, std::vector<ContigEdge> >& full_contig_graph,
                         const std::vector<ContigInfo>& candidate_contigs)
{
  std::unordered_map<std::string, std::vector<ContigEdge> > candidate_subgraph;
  
  // Build a set of candidate contig names for O(1) lookup
  std::set<std::string> candidate_contig_names;
  for(const auto& info : candidate_contigs) {
    candidate_contig_names.insert(info.contig_name);
  }
  
  // Iterate through each candidate contig and extract its edges
  for(const auto& candidate : candidate_contigs) {
    const std::string& candidate_name = candidate.contig_name;
    
    // Initialize empty edge list for this candidate (ensures all candidates appear as nodes)
    std::vector<ContigEdge> candidate_edges;
    
    // Look up this candidate's outgoing edges in the full graph
    auto it = full_contig_graph.find(candidate_name);
    if(it != full_contig_graph.end()) {
      // Filter edges: only keep edges where the sink (destination) is also a candidate
      for(const auto& edge : it->second) {
        // Only include edge if sink contig is also a candidate
        if(candidate_contig_names.find(edge.sink_contig_name) != candidate_contig_names.end()) {
          candidate_edges.push_back(edge);
        }
        // If sink is not a candidate, the edge is removed (not added to candidate_edges)
        // This handles the requirement: "remove edges from non-candidate sources to candidates"
      }
    }
    // If candidate has no outgoing edges in the full graph, candidate_edges remains empty
    
    // Add candidate to subgraph (even if it has no outgoing edges, it's still a node)
    // This ensures all candidates appear in the subgraph
    candidate_subgraph[candidate_name] = candidate_edges;
  }
  
  // Debug output: count total edges in subgraph
  std::size_t total_edges = 0;
  for(const auto& entry : candidate_subgraph) {
    total_edges += entry.second.size();
  }
    
  return candidate_subgraph;
}

/**
 * Find all paths from source nodes to sink nodes in the candidate subgraph.
 * 
 * This function:
 *   1. Identifies source nodes (nodes with no incoming edges from candidates)
 *   2. Identifies sink nodes (nodes with no outgoing edges)
 *   3. Enumerates all paths from each source to each sink using DFS
 * 
 * A path is represented as a vector of contig names in order from source to sink.
 * Cycles are avoided by tracking visited nodes during DFS.
 * 
 * @param candidate_subgraph The candidate-only subgraph adjacency list
 * @param num_source_nodes Output parameter: number of source nodes found
 * @param num_sink_nodes Output parameter: number of sink nodes found
 * @return Vector of paths, where each path is a vector of contig names (strings)
 */
std::vector<std::vector<std::string> >
build_candidate_paths(const std::unordered_map<std::string, std::vector<ContigEdge> >& candidate_subgraph,
                      std::size_t& num_source_nodes,
                      std::size_t& num_sink_nodes)
{
  std::vector<std::vector<std::string> > all_paths;
  num_source_nodes = 0;
  num_sink_nodes = 0;
  
  if(candidate_subgraph.empty()) {
    std::cout << "  No candidate contigs in subgraph, no paths to enumerate\n";
    return all_paths;
  }
  
  // Step 1: Identify source nodes (nodes with no incoming edges from candidates)
  // Build a set of all nodes that appear as sinks in edges
  std::set<std::string> nodes_with_incoming_edges;
  for(const auto& entry : candidate_subgraph) {
    for(const auto& edge : entry.second) {
      nodes_with_incoming_edges.insert(edge.sink_contig_name);
    }
  }
  
  // Source nodes are those in the subgraph but not appearing as sinks
  std::vector<std::string> source_nodes;
  for(const auto& entry : candidate_subgraph) {
    if(nodes_with_incoming_edges.find(entry.first) == nodes_with_incoming_edges.end()) {
      source_nodes.push_back(entry.first);
    }
  }
  
  // Step 2: Identify sink nodes (nodes with no outgoing edges)
  std::vector<std::string> sink_nodes;
  for(const auto& entry : candidate_subgraph) {
    if(entry.second.empty()) {
      sink_nodes.push_back(entry.first);
    }
  }
  
  // Handle case where all nodes have outgoing edges (no explicit sinks)
  // In this case, nodes with no incoming edges could be considered sources
  // and we might need to consider all nodes as potential sinks
  
  num_source_nodes = source_nodes.size();
  num_sink_nodes = sink_nodes.size();
  
  std::cout << "  Found " << source_nodes.size() << " source node(s): ";
  for(size_t i = 0; i < source_nodes.size(); ++i) {
    std::cout << source_nodes[i];
    if(i < source_nodes.size() - 1) std::cout << ", ";
  }
  std::cout << "\n";
  
  std::cout << "  Found " << sink_nodes.size() << " sink node(s): ";
  for(size_t i = 0; i < sink_nodes.size(); ++i) {
    std::cout << sink_nodes[i];
    if(i < sink_nodes.size() - 1) std::cout << ", ";
  }
  std::cout << "\n";
  
  // If no sources found, all nodes are part of cycles - treat all as potential sources
  if(source_nodes.empty()) {
    std::cout << "  Warning: No source nodes found (graph may contain cycles). Treating all nodes as potential sources.\n";
    for(const auto& entry : candidate_subgraph) {
      source_nodes.push_back(entry.first);
    }
  }
  
  // If no sinks found, all nodes have outgoing edges - treat all as potential sinks
  if(sink_nodes.empty()) {
    std::cout << "  Warning: No sink nodes found. Treating all nodes as potential sinks.\n";
    for(const auto& entry : candidate_subgraph) {
      sink_nodes.push_back(entry.first);
    }
  }
  
  // Step 3: Enumerate all paths from each source to each sink using DFS
  // Helper function for DFS path enumeration (defined as lambda with std::function for recursion)
  std::function<void(const std::string&, const std::string&, std::vector<std::string>&, std::set<std::string>&)> dfs_enumerate_paths;
  dfs_enumerate_paths = [&](const std::string& current_node,
                            const std::string& target_sink,
                            std::vector<std::string>& current_path,
                            std::set<std::string>& visited) {
    // Add current node to path and mark as visited
    current_path.push_back(current_node);
    visited.insert(current_node);
    
    // If we've reached the target sink, save this path
    if(current_node == target_sink) {
      all_paths.push_back(current_path);
    } else {
      // Continue DFS to neighbors
      auto it = candidate_subgraph.find(current_node);
      if(it != candidate_subgraph.end()) {
        for(const auto& edge : it->second) {
          // Only visit nodes that haven't been visited (avoid cycles)
          if(visited.find(edge.sink_contig_name) == visited.end()) {
            dfs_enumerate_paths(edge.sink_contig_name, target_sink, current_path, visited);
          }
        }
      }
    }
    
    // Backtrack: remove current node from path and visited set
    current_path.pop_back();
    visited.erase(current_node);
  };
  
  // Enumerate paths from each source to each sink
  for(const auto& source : source_nodes) {
    for(const auto& sink : sink_nodes) {
      if(source == sink) {
        // Single-node path
        all_paths.push_back({source});
      } else {
        std::vector<std::string> current_path;
        std::set<std::string> visited;
        dfs_enumerate_paths(source, sink, current_path, visited);
      }
    }
  }
  
  std::cout << "  Constructed " << all_paths.size() << " total path(s)\n";
  
  // Print all paths
  for(size_t i = 0; i < all_paths.size(); ++i) {
    std::cout << "    Path " << (i + 1) << ": ";
    for(size_t j = 0; j < all_paths[i].size(); ++j) {
      std::cout << all_paths[i][j];
      if(j < all_paths[i].size() - 1) std::cout << " -> ";
    }
    std::cout << "\n";
  }
  
  return all_paths;
}

/**
 * Concatenate sequences along a path, accounting for overlaps between adjacent contigs.
 * 
 * For a path like A -> B -> C, if A overlaps with B by 21bp, we concatenate:
 *   A[0:len(A)-21] + B[0:len(B)-overlap(B,C)] + C
 * 
 * This avoids duplicating the overlapping sequence during concatenation.
 * 
 * @param path Ordered list of contig names in the path
 * @param candidate_subgraph Subgraph containing edge information (overlaps)
 * @param sequences Map of contig name to nucleotide sequence
 * @return Concatenated sequence without duplicate overlaps
 */
std::string
concatenate_path_sequences(const std::vector<std::string>& path,
                          const std::unordered_map<std::string, std::vector<ContigEdge> >& candidate_subgraph,
                          const std::unordered_map<std::string, std::string>& sequences)
{
  if(path.empty()) {
    return "";
  }
  
  if(path.size() == 1) {
    // Single contig path - return full sequence
    auto it = sequences.find(path[0]);
    if(it != sequences.end()) {
      return it->second;
    }
    return "";
  }
  
  std::string concatenated;
  
  // Process each contig in the path
  for(size_t i = 0; i < path.size(); ++i) {
    const std::string& contig_name = path[i];
    auto seq_it = sequences.find(contig_name);
    if(seq_it == sequences.end()) {
      std::cerr << "Warning: contig " << contig_name << " not found in sequences\n";
      continue;
    }
    
    const std::string& full_seq = seq_it->second;
    
    if(i == 0) {
      // First contig: append full sequence (overlap will be removed when adding next)
      concatenated = full_seq;
    } else {
      // Subsequent contigs: find overlap with previous contig
      const std::string& prev_contig = path[i - 1];
      auto subgraph_it = candidate_subgraph.find(prev_contig);
      
      std::size_t overlap = 0;
      if(subgraph_it != candidate_subgraph.end()) {
        // Find the edge from previous contig to current contig
        for(const auto& edge : subgraph_it->second) {
          if(edge.sink_contig_name == contig_name) {
            overlap = edge.overlap_length;
            break;
          }
        }
      }
      
      // Remove overlap from previous contig and append current contig
      if(overlap > 0 && overlap < concatenated.size()) {
        concatenated = concatenated.substr(0, concatenated.size() - overlap);
      }
      
      // Append the current contig sequence
      concatenated += full_seq;
    }
  }
  
  return concatenated;
}

/**
 * Clip a concatenated path sequence based on minimizer offsets.
 * 
 * Clipping strategy:
 *   - First contig: clip from min_offset (with 50bp flank) to end of first contig
 *   - Middle contigs: keep entirely
 *   - Last contig: clip from start to max_offset (with 50bp flank)
 * 
 * To clip properly, we need to know where each contig starts/ends in the concatenated sequence.
 * 
 * @param concatenated_seq The concatenated sequence from concatenate_path_sequences
 * @param path Ordered list of contig names in the path
 * @param candidate_subgraph Subgraph containing edge information (overlaps)
 * @param sequences Map of contig name to nucleotide sequence
 * @param contig_info Map of contig name to ContigInfo (contains min_offset, max_offset)
 * @param keep_full_threshold Fraction threshold for keeping full sequence (default: 0.90)
 * @return Pair of (clipped_sequence, output_length)
 */
std::pair<std::string, std::size_t>
clip_path_sequence(const std::string& concatenated_seq,
                  const std::vector<std::string>& path,
                  const std::unordered_map<std::string, std::vector<ContigEdge> >& candidate_subgraph,
                  const std::unordered_map<std::string, std::string>& sequences,
                  const std::unordered_map<std::string, ContigInfo>& contig_info,
                  double keep_full_threshold = 0.90)
{
  if(path.empty() || concatenated_seq.empty()) {
    return std::make_pair("", 0);
  }
  
  // Build cumulative offsets for each contig in the concatenated sequence
  std::vector<std::pair<std::size_t, std::size_t> > contig_ranges; // (start_pos, end_pos) in concatenated_seq
  std::size_t current_pos = 0;
  
  for(size_t i = 0; i < path.size(); ++i) {
    const std::string& contig_name = path[i];
    auto seq_it = sequences.find(contig_name);
    if(seq_it == sequences.end()) {
      continue;
    }
    
    std::size_t contig_len = seq_it->second.size();
    std::size_t overlap = 0;
    
    // Calculate overlap with previous contig
    if(i > 0) {
      const std::string& prev_contig = path[i - 1];
      auto subgraph_it = candidate_subgraph.find(prev_contig);
      if(subgraph_it != candidate_subgraph.end()) {
        for(const auto& edge : subgraph_it->second) {
          if(edge.sink_contig_name == contig_name) {
            overlap = edge.overlap_length;
            break;
          }
        }
      }
    }
    
    // Adjust current position: remove overlap from previous contig
    if(i > 0 && overlap > 0) {
      current_pos -= overlap;
    }
    
    std::size_t contig_start = current_pos;
    std::size_t contig_end = current_pos + contig_len;
    
    contig_ranges.push_back(std::make_pair(contig_start, contig_end));
    current_pos = contig_end;
  }
  
  // Determine clipping positions
  const std::size_t flank_size = 50;
  std::size_t clip_start = 0;
  std::size_t clip_end = concatenated_seq.size();
  
  // Clip first contig: from min_offset (with flank) to end of first contig
  if(!path.empty() && !contig_ranges.empty()) {
    auto info_it = contig_info.find(path[0]);
    if(info_it != contig_info.end()) {
      const ContigInfo& first_info = info_it->second;
      std::size_t first_start = contig_ranges[0].first;
      std::size_t first_end = contig_ranges[0].second;
      std::size_t first_len = first_end - first_start;
      
      std::size_t first_clip_start = (first_info.min_offset >= flank_size) ?
                                     (first_info.min_offset - flank_size) : 0;
      
      // Check if we should keep full first contig
      std::size_t first_clipped_len = first_len - first_clip_start;
      if(static_cast<double>(first_clipped_len) / static_cast<double>(first_len) >= keep_full_threshold) {
        // Keep full first contig
        clip_start = first_start;
      } else {
        // Clip first contig
        clip_start = first_start + first_clip_start;
      }
    }
  }
  
  // Clip last contig: from start of last contig to max_offset (with flank)
  if(path.size() > 1 && contig_ranges.size() > 1) {
    auto info_it = contig_info.find(path.back());
    if(info_it != contig_info.end()) {
      const ContigInfo& last_info = info_it->second;
      std::size_t last_start = contig_ranges.back().first;
      std::size_t last_end = contig_ranges.back().second;
      std::size_t last_len = last_end - last_start;
      
      std::size_t last_clip_end = last_info.max_offset + flank_size + 1;
      if(last_clip_end > last_len) {
        last_clip_end = last_len;
      }
      
      // Check if we should keep full last contig
      if(static_cast<double>(last_clip_end) / static_cast<double>(last_len) >= keep_full_threshold) {
        // Keep full last contig
        clip_end = last_end;
      } else {
        // Clip last contig
        clip_end = last_start + last_clip_end;
      }
    }
  } else if(path.size() == 1) {
    // Single contig path: clip from min_offset to max_offset
    auto info_it = contig_info.find(path[0]);
    if(info_it != contig_info.end()) {
      const ContigInfo& info = info_it->second;
      std::size_t contig_len = contig_ranges[0].second - contig_ranges[0].first;
      
      std::size_t chunk_start = (info.min_offset >= flank_size) ?
                                (info.min_offset - flank_size) : 0;
      std::size_t chunk_end = info.max_offset + flank_size + 1;
      if(chunk_end > contig_len) {
        chunk_end = contig_len;
      }
      
      std::size_t clipped_len = chunk_end - chunk_start;
      if(static_cast<double>(clipped_len) / static_cast<double>(contig_len) >= keep_full_threshold) {
        clip_start = contig_ranges[0].first;
        clip_end = contig_ranges[0].second;
      } else {
        clip_start = contig_ranges[0].first + chunk_start;
        clip_end = contig_ranges[0].first + chunk_end;
      }
    }
  }
  
  // Ensure valid range
  if(clip_start >= clip_end || clip_start >= concatenated_seq.size()) {
    clip_start = 0;
  }
  if(clip_end > concatenated_seq.size()) {
    clip_end = concatenated_seq.size();
  }
  
  // Extract clipped sequence
  std::string clipped_seq = concatenated_seq.substr(clip_start, clip_end - clip_start);
  std::size_t output_len = clipped_seq.size();
  
  return std::make_pair(clipped_seq, output_len);
}

/**
 * Build PathSequence objects from candidate paths.
 * 
 * This function:
 *   1. Concatenates sequences along each path accounting for overlaps
 *   2. Clips the concatenated sequences based on minimizer offsets
 *   3. Creates PathSequence objects with appropriate labels
 * 
 * @param candidate_paths Vector of paths (each path is a vector of contig names)
 * @param candidate_subgraph Subgraph containing edge information (overlaps)
 * @param sequences Map of contig name to nucleotide sequence
 * @param contig_info Map of contig name to ContigInfo (contains min_offset, max_offset)
 * @param chunk_id Chunk identifier for labeling
 * @param keep_full_threshold Fraction threshold for keeping full sequence (default: 0.90)
 * @return Vector of PathSequence objects ready for output
 */
std::vector<PathSequence>
build_path_sequences(const std::vector<std::vector<std::string> >& candidate_paths,
                     const std::unordered_map<std::string, std::vector<ContigEdge> >& candidate_subgraph,
                     const std::unordered_map<std::string, std::string>& sequences,
                     const std::unordered_map<std::string, ContigInfo>& contig_info,
                     const std::string& chunk_id,
                     double keep_full_threshold = 0.90)
{
  std::vector<PathSequence> path_sequences;
  
  for(size_t path_idx = 0; path_idx < candidate_paths.size(); ++path_idx) {
    const auto& path = candidate_paths[path_idx];
    
    if(path.empty()) {
      continue;
    }
    
    // Step 1: Concatenate sequences accounting for overlaps
    std::string concatenated = concatenate_path_sequences(path, candidate_subgraph, sequences);
    
    if(concatenated.empty()) {
      std::cerr << "Warning: Path " << (path_idx + 1) << " produced empty concatenated sequence\n";
      continue;
    }
    
    std::size_t original_length = concatenated.size();
    
    // Step 2: Clip the concatenated sequence
    auto clipped_result = clip_path_sequence(concatenated, path, candidate_subgraph, 
                                             sequences, contig_info, keep_full_threshold);
    std::string clipped_seq = clipped_result.first;
    std::size_t output_length = clipped_result.second;
    
    // Step 3: Calculate original and clipped contig lengths
    std::vector<std::size_t> orig_lengths;
    std::vector<std::size_t> clipped_lengths;
    
    const std::size_t flank_size = 50;
    for(size_t i = 0; i < path.size(); ++i) {
      const std::string& contig_name = path[i];
      auto seq_it = sequences.find(contig_name);
      auto info_it = contig_info.find(contig_name);
      
      if(seq_it == sequences.end() || info_it == contig_info.end()) {
        orig_lengths.push_back(0);
        clipped_lengths.push_back(0);
        continue;
      }
      
      std::size_t orig_len = seq_it->second.size();
      orig_lengths.push_back(orig_len);
      
      // Calculate clipped length
      const ContigInfo& info = info_it->second;
      std::size_t clipped_len = orig_len;
      
      if(i == 0) {
        // First contig: clip from min_offset
        std::size_t clip_start = (info.min_offset >= flank_size) ?
                                 (info.min_offset - flank_size) : 0;
        std::size_t clipped_region_len = orig_len - clip_start;
        if(static_cast<double>(clipped_region_len) / static_cast<double>(orig_len) >= keep_full_threshold) {
          clipped_len = orig_len;  // Keep full
        } else {
          clipped_len = clipped_region_len;
        }
      } else if(i == path.size() - 1) {
        // Last contig: clip to max_offset
        std::size_t clip_end = info.max_offset + flank_size + 1;
        if(clip_end > orig_len) {
          clip_end = orig_len;
        }
        if(static_cast<double>(clip_end) / static_cast<double>(orig_len) >= keep_full_threshold) {
          clipped_len = orig_len;  // Keep full
        } else {
          clipped_len = clip_end;
        }
      } else {
        // Middle contigs: keep full
        clipped_len = orig_len;
      }
      
      clipped_lengths.push_back(clipped_len);
    }
    
    // Step 4: Get min_offset of start node and max_offset of sink node
    std::size_t min_offset_start = 0;
    std::size_t max_offset_sink = 0;
    
    if(!path.empty()) {
      auto start_it = contig_info.find(path[0]);
      if(start_it != contig_info.end()) {
        min_offset_start = start_it->second.min_offset;
      }
      
      auto sink_it = contig_info.find(path.back());
      if(sink_it != contig_info.end()) {
        max_offset_sink = sink_it->second.max_offset;
      }
    }
    
    // Step 5: Build label for the path
    // Join contig names with underscores (_) to avoid ambiguity with hyphens in contig IDs
    std::string label = "C" + chunk_id + "#Path" + std::to_string(path_idx + 1) + ":";
    for(size_t i = 0; i < path.size(); ++i) {
      label += path[i];
      if(i < path.size() - 1) {
        label += "_";
      }
    }
    
    // Create PathSequence object
    PathSequence path_seq{
      label,
      path,
      clipped_seq,
      original_length,
      output_length,
      orig_lengths,
      clipped_lengths,
      min_offset_start,
      max_offset_sink
    };
    
    path_sequences.push_back(path_seq);
  }
  
  return path_sequences;
}

// vector<ContigInfo>
// topological_sort_of_contigs(const std::vector<ContigInfo>& candidate_contigs, const std::unordered_map<std::string, vector<ContigEdge> >& contig_graph, std::size_t max_distance = 2)
// {
//   std::vector<std::string> output;
//   std::unordered_map<std::string, bool> visited;
//   std::stack<ContigInfo> stack;
//   for (const auto& info : candidate_contigs) {
//     if (visited.find(info.contig_name) == visited.end() || visited[info.contig_name] == false) {
//       visited[info.contig_name] = true;
      
//     }
//   }
//   set<std::string> reachable_candidate_contig_names;
// }

// bool 
// topological_sort_helper(const std::vector<ContigInfo>& candidate_contigs, const std::unordered_map<std::string, vector<ContigEdge> >& contig_graph, std::size_t max_distance = 2) {

// }



/**
 * Build labeled contig records for final output from filtered candidate contigs.
 * 
 * This function takes candidate contigs (after both minimizer threshold and
 * distance-based filtering) and creates LabeledContig records with sequences
 * and appropriate labels. It handles:
 *   1. Sequence retrieval from GFA data
 *   2. Sequence chunking (if needed) between min_offset and max_offset with flanking
 *   3. Label generation with optional offset encoding
 * 
 * Label format:
 *   - If keeping full contig: C<chunk_id>#<contig_name>
 *     Example: C7#0-1-0-0-P0
 *   - If clipped: C<chunk_id>#<contig_name>:<min_offset>-<max_offset>
 *     Example: C7#0-1-0-0-P0:278629-1137618
 * 
 * Chunking logic:
 *   - Extracts sequence from (min_offset - 50bp) to (max_offset + 50bp)
 *   - If clipped length >= keep_full_threshold * original_length, keeps full sequence
 *   - This avoids unnecessary fragmentation when minimizers cover most of the contig
 * 
 * @param candidate_contigs Vector of candidate contigs (after all filtering steps)
 * @param sequences Map of contig name to nucleotide sequence from GFA file
 * @param chunk_id Chunk identifier for labeling (e.g., "0", "1")
 * @param keep_full_threshold Fraction (0.0-1.0, default: 0.90) for keeping full contig
 *                            If clipped length >= keep_full_threshold * original_length,
 *                            keep full contig instead of clipping
 * @return Vector of LabeledContig records ready for FASTA/metadata output
 */
std::vector<LabeledContig>
build_labeled_contig_records(const std::vector<ContigInfo>& candidate_contigs,
                              const std::unordered_map<std::string, std::string>& sequences,
                              std::string chunk_id,
                              double keep_full_threshold = 0.90)
{
  std::vector<LabeledContig> output;
  
  // Process each candidate contig to create labeled output records
  for(const auto& info : candidate_contigs)
  {
    // Look up the nucleotide sequence in the GFA data
    auto iter = sequences.find(info.contig_name);
    if(iter == sequences.end())
    {
      // Warning: contig was in minimizer TSV but not found in GFA
      // This can happen if the files are out of sync
      std::cerr << "Warning: contig " << info.contig_name << " not found in GFA\n";
      continue;
    }
    const std::string& full_seq = iter->second;
    std::size_t original_len = full_seq.size();
    
    // Calculate the clipped region: between min_offset and max_offset with 50bp flanking on both sides
    // Calculate start position: min_offset - 50, but don't go below 0
    const std::size_t flank_size = 50;
    std::size_t chunk_start = (info.min_offset >= flank_size) ? 
                              (info.min_offset - flank_size) : 0;
    
    // Calculate end position: max_offset + 50 + 1 (to include the max_offset position)
    // Add 1 because offsets are 0-indexed and we want to include position max_offset
    std::size_t chunk_end = info.max_offset + flank_size + 1;
    if(chunk_end > original_len)
    {
      chunk_end = original_len;  // Don't exceed sequence length
    }
    
    // Calculate what the clipped length would be
    std::size_t clipped_len = chunk_end - chunk_start;
    
    // If the clipped length is close to the original length (>= keep_full_threshold), keep the full contig
    // This avoids unnecessary clipping when minimizers cover most of the contig
    // keep_full_threshold is configurable (default: 0.90 = 90%)
    bool keep_full_contig = (original_len > 0) && 
                            (static_cast<double>(clipped_len) / static_cast<double>(original_len) >= keep_full_threshold);
    
    std::string output_seq;
    std::size_t output_len;
    std::string label;
    
    if(keep_full_contig)
    {
      // Keep the full contig sequence without clipping
      output_seq = full_seq;
      output_len = original_len;
      // Label format without offsets when keeping full sequence
      // Example: C7#0-1-0-0-P0 (full length)
      label = "C" + chunk_id + "#" + info.contig_name;
    }
    else
    {
      // Extract the chunked sequence (between min_offset-flank and max_offset+flank)
      output_seq = full_seq.substr(chunk_start, chunk_end - chunk_start);
      output_len = output_seq.size();
      // Label format with offsets when clipping: C<chunk_id>#<contig_name>:<min_offset>-<max_offset>
      // Example: C7#0-1-0-0-P0:278629-1137618
      // This encodes the original offsets so users know what slice of the contig is being used
      label = "C" + chunk_id + "#" + info.contig_name + ":" + 
              std::to_string(info.min_offset) + "-" + 
              std::to_string(info.max_offset);
    }
    
    // Build the labeled contig record
    // original_len: full sequence length before any clipping
    // output_len: length of output sequence (full length if kept, or clipped length otherwise)
    LabeledContig contig{
      label,
      info.contig_name,
      output_seq,
      original_len,
      output_len
    };
    
    output.push_back(std::move(contig));
  }
  return output;
}

//------------------------------------------------------------------------------
// Output Functions
//------------------------------------------------------------------------------

/**
 * Write metadata TSV file with information about selected contigs.
 * 
 * Output format (tab-separated):
 *   label <tab> contig_name <tab> original_length <tab> output_length
 * 
 * The label column contains the new format (C<chunk>#<contig>), while
 * contig_name contains the original identifier from the GFA file.
 * original_length is the full contig length, while output_length is the
 * chunked length (up to max_offset where minimizers were found).
 * 
 * @param contigs Vector of labeled contig records to write
 * @param path Output file path
 * @throws std::runtime_error if file cannot be opened for writing
 */
void write_metadata(const std::vector<LabeledContig>& contigs, const std::string& path)
{
  std::ofstream out(path);
  if(!out)
  {
    throw std::runtime_error("Could not open " + path + " for writing");
  }
  
  // Write TSV header
  out << "label\tcontig_name\toriginal_length\toutput_length\n";
  
  // Write one line per contig
  for(const auto& c : contigs)
  {
    out << c.label << '\t'
        << c.contig_name << '\t'
        << c.original_len << '\t'
        << c.output_len << '\n';
  }
}


/**
 * Write FASTA file with path sequences.
 * 
 * Overloaded version for PathSequence objects. Writes concatenated and clipped
 * path sequences to FASTA format.
 * 
 * @param path_sequences Vector of PathSequence objects to write
 * @param path Output FASTA file path
 * @throws std::runtime_error if file cannot be opened for writing
 */
void write_fasta(const std::vector<PathSequence>& path_sequences, const std::string& path)
{
  std::ofstream out(path);
  if(!out)
  {
    throw std::runtime_error("Could not open " + path + " for writing");
  }
  
  // FASTA line width (standard convention: 60 characters)
  const std::size_t width = 60;
  
  // Write each path as a FASTA entry
  for(const auto& path_seq : path_sequences)
  {
    // Write header line
    out << ">" << path_seq.label << "\n";
    
    // Write sequence, wrapped at 'width' characters per line
    for(std::size_t i = 0; i < path_seq.sequence.size(); i += width)
    {
      out << path_seq.sequence.substr(i, width) << "\n";
    }
  }
}

/**
 * Write metadata TSV file for path sequences.
 * 
 * Overloaded version for PathSequence objects. Output format:
 *   label <tab> path_contigs <tab> original_contig_lengths <tab> clipped_contig_lengths
 *   <tab> min_offset_in_start_node <tab> max_offset_in_sink_node
 * 
 * @param path_sequences Vector of PathSequence objects to write
 * @param path Output TSV file path
 * @throws std::runtime_error if file cannot be opened for writing
 */
void write_metadata(const std::vector<PathSequence>& path_sequences, const std::string& path)
{
  std::ofstream out(path);
  if(!out)
  {
    throw std::runtime_error("Could not open " + path + " for writing");
  }
  
  // Write TSV header
  out << "label\tpath_contigs\toriginal_contig_lengths\tclipped_contig_lengths\tmin_offset_in_start_node\tmax_offset_in_sink_node\n";
  
  // Write one line per path
  for(const auto& path_seq : path_sequences)
  {
    // Build path contigs string (comma-separated)
    std::string path_contigs_str;
    for(size_t i = 0; i < path_seq.path_contigs.size(); ++i) {
      path_contigs_str += path_seq.path_contigs[i];
      if(i < path_seq.path_contigs.size() - 1) {
        path_contigs_str += ",";
      }
    }
    
    // Build original contig lengths string (comma-separated)
    std::string orig_lengths_str;
    for(size_t i = 0; i < path_seq.original_contig_lengths.size(); ++i) {
      orig_lengths_str += std::to_string(path_seq.original_contig_lengths[i]);
      if(i < path_seq.original_contig_lengths.size() - 1) {
        orig_lengths_str += ",";
      }
    }
    
    // Build clipped contig lengths string (comma-separated)
    std::string clipped_lengths_str;
    for(size_t i = 0; i < path_seq.clipped_contig_lengths.size(); ++i) {
      clipped_lengths_str += std::to_string(path_seq.clipped_contig_lengths[i]);
      if(i < path_seq.clipped_contig_lengths.size() - 1) {
        clipped_lengths_str += ",";
      }
    }
    
    out << path_seq.label << '\t'
        << path_contigs_str << '\t'
        << orig_lengths_str << '\t'
        << clipped_lengths_str << '\t'
        << path_seq.min_offset_in_start_node << '\t'
        << path_seq.max_offset_in_sink_node << '\n';
  }
}

/**
 * Write TSV file with distribution data for downstream plotting.
 * 
 * Creates a TSV file containing all contigs with their minimizer hit counts,
 * selection status, and metadata. This can be used by external tools (R, Python,
 * etc.) to generate distribution plots.
 * 
 * Output format (tab-separated):
 *   contig_name <tab> label <tab> minimizer_hits <tab> selected <tab> 
 *   original_length <tab> output_length <tab> min_offset <tab> max_offset
 * 
 * The 'selected' column indicates whether the contig passed the percentile threshold.
 * 'original_length' is the full contig length, 'output_length' is the chunked length
 * (up to max_offset). Offsets indicate where minimizers were found on the contig.
 * 
 * @param minimizer_counts Vector of all minimizer hit counts (for statistics, currently unused but kept for future use)
 * @param all_contigs Map of all contig metadata
 * @param selected_contigs Vector of contigs that passed the threshold filter
 * @param threshold The calculated threshold value (percentile or IQR-based)
 * @param method_value The method parameter (percentile value or IQR multiplier)
 * @param plot_data_path Path where TSV data file should be written
 * @param use_iqr Whether IQR method was used (true) or percentile method (false)
 * @throws std::runtime_error if file cannot be opened for writing
 */
void write_distribution_data(const std::vector<std::size_t>& minimizer_counts,
                            const std::unordered_map<std::string, ContigInfo>& all_contigs,
                            const std::vector<LabeledContig>& selected_contigs,
                            std::size_t threshold,
                            double method_value,
                            const std::string& plot_data_path,
                            bool use_iqr = false)
{
  // Open TSV file for writing
  std::ofstream data_out(plot_data_path);
  if(!data_out)
  {
    throw std::runtime_error("Could not open " + plot_data_path + " for writing");
  }
  
  // Build a hash set of selected contig names for O(1) lookup
  std::unordered_set<std::string> selected_contig_names;
  // Also build a map from contig name to LabeledContig for quick access
  std::unordered_map<std::string, const LabeledContig*> selected_map;
  for(const auto& c : selected_contigs)
  {
    selected_contig_names.insert(c.contig_name);
    selected_map[c.contig_name] = &c;
  }
  
  // Write TSV header (tab-separated)
  data_out << "contig_name\tlabel\tminimizer_hits\tselected\toriginal_length\toutput_length\tmin_offset\tmax_offset\n";
  
  // Write all contigs with their minimizer hit counts, selection status, and metadata
  for(const auto& item : all_contigs)
  {
    const ContigInfo& info = item.second;
    
    // Check if this contig was selected (passed the threshold)
    bool is_selected = selected_contig_names.count(info.contig_name) > 0;
    
    // Retrieve label and length info if contig was selected
    std::string label = "";
    std::size_t original_length = 0;
    std::size_t output_length = 0;
    
    if(is_selected)
    {
      auto map_iter = selected_map.find(info.contig_name);
      if(map_iter != selected_map.end())
      {
        label = map_iter->second->label;
        original_length = map_iter->second->original_len;
        output_length = map_iter->second->output_len;
      }
    }
    
    // Write TSV row (tab-separated):
    // contig_name, label, minimizer_hits, selected (TRUE/FALSE), 
    // original_length, output_length, min_offset, max_offset
    data_out << info.contig_name << "\t"
             << label << "\t"
             << info.unique_minimizers << "\t"
             << (is_selected ? "TRUE" : "FALSE") << "\t"
             << original_length << "\t"
             << output_length << "\t"
             << info.min_offset << "\t"
             << info.max_offset << "\n";
  }
  
  data_out.close();
}

//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------

int main(int argc, char** argv)
{
  if(argc < 11)
  {
    std::cerr << "Usage: " << argv[0]
              << " --chunk-id <id>"
              << " --contig-type <target|query>"
              << " --contig-tsv <file>"
              << " --gfa <assembly.gfa>"
              << " [--percentile <N> | --iqr <multiplier>]"
              << " --fasta-output <paths.fasta>"
              << " --metadata-output <paths.tsv>"
              << " --stats-output <stats.tsv>\n"
              << "\n"
              << "This tool extracts candidate paths from contigs that are outliers in minimizer hit distribution.\n"
              << "Contigs can be selected using either percentile or IQR (Interquartile Range) method.\n"
              << "Output includes concatenated and clipped path sequences:\n"
              << "  - Paths are built from candidate contigs connected in the assembly graph\n"
              << "  - Sequences are concatenated accounting for overlaps between contigs\n"
              << "  - First and last contigs are clipped based on minimizer offsets (with 50bp flanking)\n"
              << "  - Middle contigs are kept entirely\n"
              << "  - If clipped length is close to original (see --keep-full-threshold), keep full contig\n"
              << "\n"
             << "Label format: C<chunk>#Path<path_id>:<contig1>_<contig2>_...\n"
             << "Example: C7#Path1:0-1-0-0-P0_1-4-0-0-P1\n"
              << "\n"
              << "Selection Methods (use one of the following):\n"
              << "  --percentile <N>     Percentile threshold (0-100, default: 90)\n"
              << "                       Contigs >= Nth percentile are selected\n"
              << "  --iqr <multiplier>   IQR-based outlier detection (e.g., 1.5 or 3.0)\n"
              << "                       Contigs > Q3 + (multiplier * IQR) are selected\n"
              << "                       Use 1.5 for mild outliers, 3.0 for extreme outliers\n"
              << "\n"
              << "Optional Arguments:\n"
              << "  --keep-full-threshold <fraction>  Fraction (0.0-1.0, default: 0.90)\n"
              << "                                    If clipped length >= fraction * original length,\n"
              << "                                    keep full contig instead of clipping\n";
    return 1;
  }

  std::string chunk_id;
  std::string contig_type;
  std::string contig_tsv;
  std::string gfa_path;
  double percentile = -1.0;  // -1 indicates not set, use IQR instead
  double iqr_multiplier = -1.0;  // -1 indicates not set, use percentile instead
  double keep_full_threshold = 0.90;  // Default: 90% - keep full contig if clipped length >= 90% of original
  std::string fasta_output;
  std::string metadata_output;
  std::string stats_output;

  for(int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];
    if(arg == "--chunk-id" && i + 1 < argc)
    {
      chunk_id = argv[++i];
    }
    else if(arg == "--contig-type" && i + 1 < argc)
    {
      contig_type = argv[++i];
    }
    else if(arg == "--contig-tsv" && i + 1 < argc)
    {
      contig_tsv = argv[++i];
    }
    else if(arg == "--gfa" && i + 1 < argc)
    {
      gfa_path = argv[++i];
    }
    else if(arg == "--percentile" && i + 1 < argc)
    {
      percentile = std::stod(argv[++i]);
      if(percentile < 0.0 || percentile > 100.0)
      {
        std::cerr << "Error: percentile must be between 0 and 100\n";
        return 1;
      }
    }
    else if(arg == "--iqr" && i + 1 < argc)
    {
      iqr_multiplier = std::stod(argv[++i]);
      if(iqr_multiplier < 0.0)
      {
        std::cerr << "Error: IQR multiplier must be >= 0\n";
        return 1;
      }
    }
    else if(arg == "--fasta-output" && i + 1 < argc)
    {
      fasta_output = argv[++i];
    }
    else if(arg == "--metadata-output" && i + 1 < argc)
    {
      metadata_output = argv[++i];
    }
    else if(arg == "--stats-output" && i + 1 < argc)
    {
      stats_output = argv[++i];
    }
    else if(arg == "--keep-full-threshold" && i + 1 < argc)
    {
      keep_full_threshold = std::stod(argv[++i]);
      if(keep_full_threshold < 0.0 || keep_full_threshold > 1.0)
      {
        std::cerr << "Error: keep-full-threshold must be between 0.0 and 1.0\n";
        return 1;
      }
    }
    else
    {
      std::cerr << "Unknown or incomplete argument: " << arg << "\n";
      return 1;
    }
  }

  if(chunk_id.empty() || contig_type.empty() || contig_tsv.empty() || gfa_path.empty() || 
     fasta_output.empty() || metadata_output.empty() || stats_output.empty())
  {
    std::cerr << "Missing required arguments\n";
    return 1;
  }
  
  // Validate contig_type
  if(contig_type != "target" && contig_type != "query")
  {
    std::cerr << "Error: --contig-type must be either 'target' or 'query'\n";
    return 1;
  }
  
  // Validate that exactly one selection method is specified
  bool use_percentile = (percentile >= 0.0);
  bool use_iqr = (iqr_multiplier >= 0.0);
  
  if(!use_percentile && !use_iqr)
  {
    // Default to percentile method if neither specified
    percentile = 90.0;
    use_percentile = true;
  }
  else if(use_percentile && use_iqr)
  {
    std::cerr << "Error: Cannot use both --percentile and --iqr. Please specify only one selection method.\n";
    return 1;
  }

  try
  {
    // ========================================================================
    // Step 1: Load minimizer statistics from TSV file
    // ========================================================================
    std::cout << "================================================================" << "\n";
    std::cout << "Loading contig minimizer info from " << contig_tsv << "\n";
    std::cout << "================================================================" << "\n";
    auto contig_info = load_contig_filter(contig_tsv);
    std::cout << "  Found " << contig_info.size() << " contigs in TSV\n";
    
    if(contig_info.empty())
    {
      std::cerr << "Error: No contigs found in TSV file\n";
      return 1;
    }

    // ========================================================================
    // Step 2: Build contig adjacency graph
    // ========================================================================
    auto contig_adjacency_graph = build_contig_graph(gfa_path);
    
    // ========================================================================
    // Step 3: Calculate distribution and percentile threshold
    // ========================================================================
    // Extract all minimizer hit counts into a vector for analysis
    std::vector<std::size_t> minimizer_counts;
    minimizer_counts.reserve(contig_info.size());
    for(const auto& item : contig_info)
    {
      minimizer_counts.push_back(item.second.unique_minimizers);
    }
    
    // Sort counts for threshold calculation (requires sorted input)
    std::sort(minimizer_counts.begin(), minimizer_counts.end());
    
    // Calculate threshold using selected method
    std::size_t threshold;
    std::string method_description;
    
    if(use_iqr)
    {
      // Calculate IQR-based threshold
      threshold = calculate_iqr_threshold(minimizer_counts, iqr_multiplier);
      double q1 = calculate_percentile(minimizer_counts, 25.0);
      double q3 = calculate_percentile(minimizer_counts, 75.0);
      double iqr = static_cast<double>(q3) - static_cast<double>(q1);
      
      method_description = "IQR method (multiplier: " + std::to_string(iqr_multiplier) + ")";
      
      // Print distribution statistics
      std::cout << "\nMinimizer hit distribution:\n";
      std::cout << "  Total contigs: " << minimizer_counts.size() << "\n";
      std::cout << "  Minimum: " << minimizer_counts.front() << "\n";
      std::cout << "  Maximum: " << minimizer_counts.back() << "\n";
      std::cout << "  Q1 (25th percentile): " << static_cast<std::size_t>(q1) << "\n";
      std::cout << "  Median (50th percentile): " << calculate_percentile(minimizer_counts, 50.0) << "\n";
      std::cout << "  Q3 (75th percentile): " << static_cast<std::size_t>(q3) << "\n";
      std::cout << "  IQR: " << iqr << "\n";
      std::cout << "  Upper fence (Q3 + " << iqr_multiplier << " * IQR): " << threshold << "\n";
      std::cout << "  Using threshold: >= " << threshold << " minimizer hits\n";
    }
    else
    {
      // Calculate percentile-based threshold
      threshold = calculate_percentile(minimizer_counts, percentile);
      method_description = std::to_string(static_cast<int>(percentile)) + "th percentile";
      
      // Print distribution statistics
      std::cout << "\nMinimizer hit distribution:\n";
      std::cout << "  Total contigs: " << minimizer_counts.size() << "\n";
      std::cout << "  Minimum: " << minimizer_counts.front() << "\n";
      std::cout << "  Maximum: " << minimizer_counts.back() << "\n";
      std::cout << "  Median: " << calculate_percentile(minimizer_counts, 50.0) << "\n";
      std::cout << "  " << percentile << "th percentile: " << threshold << "\n";
      std::cout << "  Using threshold: >= " << threshold << " minimizer hits\n";
    }
    
    // ========================================================================
    // Step 4: Load nucleotide sequences from GFA file
    // ========================================================================
    std::cout << "================================================================" << "\n";
    std::cout << "Loading GFA sequences from " << gfa_path << "\n";
    std::cout << "================================================================" << "\n";
    auto sequences = load_gfa_sequences(gfa_path);
    std::cout << "  Found " << sequences.size() << " sequences in GFA\n";
    
    // ========================================================================
    // Step 5: Filter contigs and build output records
    // ========================================================================
    // Select contigs that meet or exceed the threshold
    // These are considered outliers with high minimizer hit frequency
    auto candidate_configs = filter_candidate_contig_records_by_minimizer_threshold(contig_info, threshold, keep_full_threshold);
    
    // Collect statistics for TSV output
    std::size_t num_initial_contigs = contig_info.size();
    std::size_t num_contigs_after_threshold = candidate_configs.size();
    std::vector<std::size_t> component_minimizer_totals;
    std::vector<std::size_t> component_contig_totals;
    
    auto filtered_candidate_configs = select_connected_candidate_contigs(contig_adjacency_graph, candidate_configs, 2, contig_type, component_minimizer_totals, component_contig_totals);
    
    std::size_t num_contigs_in_selected_component = filtered_candidate_configs.size();

    
    std::cout << "\n  candidate_contigs:\n";
    for(const auto& info : filtered_candidate_configs)
    {
      std::cout << "  " << info.contig_name << "\n";
    }
    
    // ========================================================================
    // Step 6: Build candidate subgraph and candidate paths
    // ========================================================================
    std::cout << "================================================================\n";
    std::cout << "BUILDING CANDIDATE PATHS\n";
    std::cout << "================================================================\n";
    
    /*
    TODO: When creating candidate subgraph, find start and sink nodes, and extract (from adjacency graph) the entire subgraph starting from start node and ending at sink node. 
    This will ensure that any incorrectly filtered out candidate contigs that are in overlap region (i.e., between start and sink nodes) are included in the subgraph.
    */
    auto candidate_subgraph = build_candidate_subgraph(contig_adjacency_graph, filtered_candidate_configs);
    // Calculate total edges in the candidate subgraph
    std::size_t total_edges = 0;
    for(const auto& entry : candidate_subgraph) {
      total_edges += entry.second.size();
    }
    std::cout << "  Built candidate subgraph with " << candidate_subgraph.size() 
              << " candidate contig nodes and " << total_edges 
              << " edges between candidate contigs\n\n";

    std::size_t num_source_nodes = 0;
    std::size_t num_sink_nodes = 0;
    auto candidate_paths = build_candidate_paths(candidate_subgraph, num_source_nodes, num_sink_nodes);
    
    std::size_t num_candidate_paths = candidate_paths.size();
    
    // ========================================================================
    // Step 7: Build path sequences (concatenate and clip)
    // ========================================================================
    std::cout << "================================================================\n";
    std::cout << "BUILDING PATH SEQUENCES (CONCATENATION AND CLIPPING)\n";
    std::cout << "================================================================\n";
    auto path_sequences = build_path_sequences(candidate_paths, candidate_subgraph, 
                                               sequences, contig_info, chunk_id, keep_full_threshold);
    
    std::cout << "  Built " << path_sequences.size() << " path sequence(s)\n";
    for(size_t i = 0; i < path_sequences.size(); ++i) {
      std::cout << "    Path " << (i + 1) << ": " << path_sequences[i].label 
                << " (original: " << path_sequences[i].original_length 
                << " bp, clipped: " << path_sequences[i].output_length << " bp)\n";
    }
    std::cout << "\n";

    // ========================================================================
    // Step 8: Write output files
    // ========================================================================
    std::cout << "================================================================\n";
    std::cout << "WRITING OUTPUT FILES (FASTA, METADATA)\n";
    std::cout << "================================================================\n";
    // Write FASTA file with path sequences
    write_fasta(path_sequences, fasta_output);
    
    // Write TSV metadata file with path information
    write_metadata(path_sequences, metadata_output);
    
    // Write statistics TSV file
    std::ofstream stats_out(stats_output);
    if(!stats_out)
    {
      throw std::runtime_error("Could not open " + stats_output + " for writing");
    }
    
    // Write TSV header
    stats_out << "num_initial_contigs\tnum_contigs_after_threshold\tcomponent_minimizer_totals\tcomponent_contig_totals\tnum_contigs_in_selected_component\tnum_source_sink_nodes\tnum_candidate_paths\n";
    
    // Build component minimizer totals string (tuple format)
    std::string component_minimizer_totals_str = "(";
    for(size_t i = 0; i < component_minimizer_totals.size(); ++i) {
      component_minimizer_totals_str += std::to_string(component_minimizer_totals[i]);
      if(i < component_minimizer_totals.size() - 1) {
        component_minimizer_totals_str += ",";
      }
    }
    component_minimizer_totals_str += ")";
    
    // Build component contig totals string (tuple format)
    std::string component_contig_totals_str = "(";
    for(size_t i = 0; i < component_contig_totals.size(); ++i) {
      component_contig_totals_str += std::to_string(component_contig_totals[i]);
      if(i < component_contig_totals.size() - 1) {
        component_contig_totals_str += ",";
      }
    }
    component_contig_totals_str += ")";
    
    // Build source/sink nodes tuple
    std::string source_sink_tuple = "(" + std::to_string(num_source_nodes) + "," + std::to_string(num_sink_nodes) + ")";
    
    // Write statistics row
    stats_out << num_initial_contigs << "\t"
              << num_contigs_after_threshold << "\t"
              << component_minimizer_totals_str << "\t"
              << component_contig_totals_str << "\t"
              << num_contigs_in_selected_component << "\t"
              << source_sink_tuple << "\t"
              << num_candidate_paths << "\n";
    
    stats_out.close();
    
    // Print summary of generated files
    std::cout << "\nWrote FASTA to " << fasta_output << "\n";
    std::cout << "Wrote metadata to " << metadata_output << "\n";
    std::cout << "Wrote statistics to " << stats_output << "\n";
  }
  catch(const std::exception& ex)
  {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }

  return 0;
}

