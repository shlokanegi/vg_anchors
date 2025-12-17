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
  
  // Debug output: print filtered candidate contigs
  std::cout << "filtered candidate_contigs by minimizer threshold. Num remaining: " << candidate_contigs.size() << std::endl;
  std::cout << "candidate_contigs: " << std::endl;
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
 * This graph is used by filter_distant_candidate_contigs to perform distance-based
 * filtering (BFS) to find candidate contigs that are within max_distance hops of
 * each other.
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
    // Only process L-lines (link/overlap lines)
    if(line[0] == 'L')
    {
      std::string source_contig_name, sink_contig_name;
      std::string link_line_symbol;  // Will be 'L'
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

/**
 * Filter candidate contigs to keep only those within max_distance hops of each other.
 * 
 * This function performs distance-based filtering using BFS (Breadth-First Search)
 * on the contig adjacency graph. The goal is to keep candidate contigs that form
 * "clusters" or connected components within the graph, filtering out isolated
 * candidates that are too far from other candidates.
 * 
 * Algorithm:
 *   1. For each candidate contig, perform BFS starting from that contig
 *   2. Explore up to max_distance hops in the graph
 *   3. If another candidate contig is found within max_distance, mark both the
 *      starting candidate and the found candidate as "reachable"
 *   4. Only candidates that are reachable from at least one other candidate
 *      (or are part of a cluster) are kept in the output
 * 
 * Special case: If there's only 1 candidate, it's automatically kept (no need
 * to check reachability to other candidates).
 * 
 * This filtering step helps remove isolated candidate contigs that may be false
 * positives or not part of the same genomic region as other candidates.
 * 
 * @param contig_graph Adjacency list representation of contig overlaps (from GFA L-lines)
 * @param candidate_contigs Vector of candidate contigs (already filtered by minimizer threshold)
 * @param max_distance Maximum number of hops allowed between candidate contigs (default: 2)
 *                    Contigs more than max_distance apart are considered too distant
 * @return Vector of candidate contigs that are within max_distance of at least one other candidate
 */
std::vector<ContigInfo>
filter_distant_candidate_contigs(const std::unordered_map<std::string, std::vector<ContigEdge> >& contig_graph,
                                 const std::vector<ContigInfo>& candidate_contigs,
                                 std::size_t max_distance = 2)
{
  std::vector<ContigInfo> output;
  
  // Build a set of candidate contig names for O(1) lookup during BFS
  std::set<std::string> candidate_contig_names;
  for(const auto& info : candidate_contigs) {
    candidate_contig_names.insert(info.contig_name);
  }
  
  // Track which candidate contigs are reachable from other candidates
  std::set<std::string> reachable_candidate_contig_names;
  
  // Special case: if only one candidate, keep it automatically
  if (candidate_contig_names.size() == 1) {
    reachable_candidate_contig_names.insert(candidate_contigs[0].contig_name);
  }
  else {
    // For each candidate contig, perform BFS to find other candidates within max_distance
    for (const auto& curr_contig: candidate_contigs) {
      std::string start_contig_name = curr_contig.contig_name;  // Starting point for this BFS
      
      // BFS queue: stores (contig_name, distance_from_start)
      std::queue<std::pair<std::string, std::size_t> > bfs_queue;
      bfs_queue.push(std::make_pair(start_contig_name, 0));
      
      // BFS traversal: explore graph up to max_distance hops
      while (!bfs_queue.empty()) {
        std::string curr_iter_contig_name = bfs_queue.front().first;
        std::size_t curr_distance = bfs_queue.front().second;
        bfs_queue.pop();
        
        // Stop exploring if we've exceeded max_distance
        if (curr_distance >= max_distance) {
          continue;
        }
        
        // Check if this contig has outgoing edges in the graph
        auto it = contig_graph.find(curr_iter_contig_name);
        if (it != contig_graph.end()) {
          // Explore all neighbors (sink contigs) of current contig
          for (const auto& edge : it->second) {
            // Check if the neighbor is itself a candidate contig
            if (candidate_contig_names.find(edge.sink_contig_name) != candidate_contig_names.end()) {
              // Found another candidate within reach! Mark both as reachable:
              // - The found candidate (edge.sink_contig_name)
              // - The starting candidate (start_contig_name)
              reachable_candidate_contig_names.insert(edge.sink_contig_name);
              reachable_candidate_contig_names.insert(start_contig_name);
              // Note: We don't continue BFS from found candidates here to avoid
              // redundant searches (they will be processed when we start BFS from them)
            }
            else {
              // Not a candidate, but continue BFS to potentially find candidates further away
              bfs_queue.push(std::make_pair(edge.sink_contig_name, curr_distance + 1));
            }
          }
        }
      }
    }
  }
  
  // Filter output to only include reachable candidates
  for (const auto& info : candidate_contigs) {
    if (reachable_candidate_contig_names.find(info.contig_name) != reachable_candidate_contig_names.end()) {
      output.push_back(info);
    }
  }
  
  // Debug output
  std::cout << "filtered distant candidate_contigs. Num remaining: " << output.size() << std::endl;

  return output;
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
 * Write FASTA file with selected contig sequences.
 * 
 * Standard FASTA format:
 *   >sequence_label
 *   SEQUENCE_DATA (wrapped at 60 characters per line)
 * 
 * Uses the new label format (C<chunk>#<contig>:<min_offset>-<max_offset>)
 * as the sequence header, encoding the offsets to show what slice of the
 * original contig is being used. Sequences are wrapped at 60 characters per
 * line as per FASTA convention.
 * 
 * @param contigs Vector of labeled contig records to write
 * @param path Output FASTA file path
 * @throws std::runtime_error if file cannot be opened for writing
 */
void write_fasta(const std::vector<LabeledContig>& contigs, const std::string& path)
{
  std::ofstream out(path);
  if(!out)
  {
    throw std::runtime_error("Could not open " + path + " for writing");
  }
  
  // FASTA line width (standard convention: 60 characters)
  const std::size_t width = 60;
  
  // Write each contig as a FASTA entry
  for(const auto& c : contigs)
  {
    // Write header line
    out << ">" << c.label << "\n";
    
    // Write sequence, wrapped at 'width' characters per line
    for(std::size_t i = 0; i < c.sequence.size(); i += width)
    {
      out << c.sequence.substr(i, width) << "\n";
    }
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
              << " --contig-tsv <file>"
              << " --gfa <assembly.gfa>"
              << " [--percentile <N> | --iqr <multiplier>]"
              << " --fasta-output <contigs.fasta>"
              << " --metadata-output <contigs.tsv>"
              << " --distribution-tsv <distribution.tsv>\n"
              << "\n"
              << "This tool extracts individual contigs that are outliers in minimizer hit distribution.\n"
              << "Contigs can be selected using either percentile or IQR (Interquartile Range) method.\n"
              << "Output includes chunked contig sequences (between min/max minimizer offsets with\n"
              << "50bp flanking on both sides) and a TSV file for downstream plotting.\n"
              << "If the clipped length is close to the original length (see --keep-full-threshold),\n"
              << "the full contig is kept instead of clipping.\n"
              << "\n"
              << "Label format: C<chunk>#<contig_name>:<min_offset>-<max_offset> (when clipped)\n"
              << "              C<chunk>#<contig_name> (when keeping full sequence)\n"
              << "Example: C7#0-1-0-0-P0:278629-1137618 (clipped) or C7#0-1-0-0-P0 (full)\n"
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
  std::string contig_tsv;
  std::string gfa_path;
  double percentile = -1.0;  // -1 indicates not set, use IQR instead
  double iqr_multiplier = -1.0;  // -1 indicates not set, use percentile instead
  double keep_full_threshold = 0.90;  // Default: 90% - keep full contig if clipped length >= 90% of original
  std::string fasta_output;
  std::string metadata_output;
  std::string distribution_tsv;

  for(int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];
    if(arg == "--chunk-id" && i + 1 < argc)
    {
      chunk_id = argv[++i];
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
    else if(arg == "--distribution-tsv" && i + 1 < argc)
    {
      distribution_tsv = argv[++i];
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

  if(chunk_id.empty() || contig_tsv.empty() || gfa_path.empty() || 
     fasta_output.empty() || metadata_output.empty() || distribution_tsv.empty())
  {
    std::cerr << "Missing required arguments\n";
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
    std::cout << "Loading contig minimizer info from " << contig_tsv << "\n";
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
    std::cout << "\nLoading GFA sequences from " << gfa_path << "\n";
    auto sequences = load_gfa_sequences(gfa_path);
    std::cout << "  Found " << sequences.size() << " sequences in GFA\n";
    
    // ========================================================================
    // Step 5: Filter contigs and build output records
    // ========================================================================
    // Select contigs that meet or exceed the threshold
    // These are considered outliers with high minimizer hit frequency
    auto candidate_configs = filter_candidate_contig_records_by_minimizer_threshold(contig_info, threshold, keep_full_threshold);
    
    auto filtered_candidate_configs = filter_distant_candidate_contigs(contig_adjacency_graph, candidate_configs, 2);
    auto records = build_labeled_contig_records(filtered_candidate_configs, sequences, chunk_id, keep_full_threshold);

    std::cout << "\nExtracted " << records.size() << " contigs (>= " << threshold 
              << " minimizer hits, " << method_description << "):\n";
    for(const auto& r : records)
    {
      std::cout << "  " << r.contig_name << "\n";
    }
    
    // ========================================================================
    // Step 6: Write output files
    // ========================================================================
    // Write FASTA file with chunked contig sequences (up to max_offset)
    write_fasta(records, fasta_output);
    
    // Write TSV metadata file with contig information
    write_metadata(records, metadata_output);
    
    // Write TSV file with distribution data for downstream plotting
    // Pass the appropriate threshold and method info
    double method_value = use_iqr ? iqr_multiplier : percentile;
    write_distribution_data(minimizer_counts, contig_info, records, 
                           threshold, method_value, distribution_tsv, use_iqr);
    
    // Print summary of generated files
    std::cout << "\nWrote FASTA to " << fasta_output << "\n";
    std::cout << "Wrote metadata to " << metadata_output << "\n";
    std::cout << "Wrote distribution TSV to " << distribution_tsv << "\n";
  }
  catch(const std::exception& ex)
  {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }

  return 0;
}

