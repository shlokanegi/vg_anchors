from assembler.gaf_reader import GafReader
#from assembler.builder import AnchorDictionary
from assembler.aligner import AlignAnchor
import assembler.parser as parser
import time
import multiprocessing
from sys import stderr
import os
from assembler.config import settings
from collections import defaultdict
from assembler.read import Read
from memory_profiler import profile as mem_profile
from line_profiler import profile as line_profile
import sys

# Make @profile available - it acts as a no-op when not using kernprof
try:
    from line_profiler import profile
except ImportError:
    # If line_profiler is not available, create a no-op decorator
    def profile(func):
        return func

# Ensure we use fork mode for true copy-on-write behavior
# (on Linux, this is the default, but we make it explicit for clarity)
try:
    multiprocessing.set_start_method('fork', force=False)
except RuntimeError:
    # Already set, which is fine
    pass

# Global object to hold shared data for worker processes
# Note: Set this before creating the multiprocessing pool to leverage fork()'s copy-on-write
shared_align_anchor = None

def nested_dd_factory():
    return defaultdict(list)

def init_worker():
    """
    Initializer for each worker process in the pool.
    On Linux with fork, the global shared_align_anchor is already available via copy-on-write.
    """
    # No need to do anything - fork() gives us access to the parent's memory
    pass


def process_gaf_chunk(gaf_chunk_lines: list[str]) -> dict:
    """
    Worker function to process a chunk of GAF lines.
    This function is executed in a separate process.
    """
    global shared_align_anchor

    # Initialize local dictionaries to store results for this chunk.
    local_anchor_reads_dict = defaultdict(nested_dd_factory)
    local_bp_matched_reads = defaultdict(list)
    local_read_ranks = defaultdict(list)
    local_path_matched_reads = defaultdict(list)
    local_reads_processed_dict = {} # {read_name: processed_line_data}

    t0 = time.time()
    # Process each line in the assigned GAF chunk.
    for line in gaf_chunk_lines:
        processed_line_data = parser.processGafLine(line)

        if shared_align_anchor.read_id_map:
            if not processed_line_data:
                continue
            read_name = processed_line_data[0]
            read_id = shared_align_anchor.read_id_map.get(read_name)
            if read_id is None:
                continue # Skip GAF entries for reads not in the FASTA file
            processed_line_data[0] = read_id

        if processed_line_data:
            if settings.OUTPUT_LOGGING_FILES:
                local_reads_processed_dict[processed_line_data[0]] = processed_line_data
            # remove mapq (index MAP_Q_ID) and div (index DIV_ID) from processed_line_data
            processed_line_data = processed_line_data[:4] + processed_line_data[6:]
        
            # Call the refactored processGafLine on the shared object
            # This is a read-only operation on shared_align_anchor
            result, current_read = shared_align_anchor.processGafLine(processed_line_data)
                        
            if settings.MIN_ANCHOR_LENGTH == 0:
                for (sentinel, i), reads in result["path_matched_reads"].items():
                    local_path_matched_reads[(sentinel, i)].extend(reads)

            else:
                for (sentinel, i), reads in result["anchor_reads"].items():
                    local_anchor_reads_dict[sentinel][i].extend(reads)
                
                for (sentinel, i), reads in result["bp_matched_reads"].items():
                    local_bp_matched_reads[(sentinel, i)].extend(reads)
                
                for (sentinel, i), read_rank in result["read_ranks"].items():
                    local_read_ranks[(sentinel, i)].extend(read_rank)
                
                if settings.OUTPUT_LOGGING_FILES:
                    for (sentinel, i), reads in result["path_matched_reads"].items():
                        local_path_matched_reads[(sentinel, i)].extend(reads)

    if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
        print(f" ..Processed {len(gaf_chunk_lines)} lines in {time.time()-t0:.2f}s", file=stderr)

    # Return the collected results from this worker.
    if settings.MIN_ANCHOR_LENGTH == 0:
        return {
            "path_matched_reads": local_path_matched_reads,
            "reads_processed": local_reads_processed_dict
        }
    else:
        return {
            "anchor_reads_dict": local_anchor_reads_dict,
            "bp_matched_reads": local_bp_matched_reads,
            "read_ranks": local_read_ranks,
            "path_matched_reads": local_path_matched_reads,
            "reads_processed": local_reads_processed_dict
        }


class Orchestrator:

    def __init__(
        self, dictionary_path: str, graph_path: str, gaf_path: str, fasta_path: str, threads: int, read_id_map: dict = None
    ):
        """
        It initiailzes the AlignAnchor object with the packedgraph path and the dictionary generated by the assembler.builder.AnchorDictionrary object.
        It initializes the GafReader object that reads the gaf file.

        Parameters
        ----------
        sentinel_to_anchor_dictionary: dictionary
            the dctionary associating sentinels and anchors
        graph_path: string
            The filepath of the packedGraph object
        gaf_path:
            The filepath of the gaf alignment file
        fasta_path: string
            The filepath of the reads fasta file
        """
        self.align_anchor = AlignAnchor(threads=int(threads), read_id_map=read_id_map)
        t0 = time.time()
        self.align_anchor.build(dictionary_path, graph_path)    # graph is loaded here once!
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(f"AlignAnchor built in {time.time()-t0:.2f}s", file=stderr)
        self.align_anchor.readFasta(fasta_path)
        self.gaf_path = gaf_path
        self.threads = int(threads)
        self.total_reads_in_gaf = 0

    def _chunk_gaf_file(self, gaf_path: str, num_chunks: int) -> list:
        """
        Reads a GAF file and splits its lines into a specified number of chunks for parallel processing.
        """
        with open(gaf_path, "r") as f:
            lines = f.readlines()
        
        self.total_reads_in_gaf = len(lines)
        if not lines:
            return []

        chunk_size = (len(lines) + num_chunks - 1) // num_chunks
        return [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)] # output: [[line1, line2, ...], [line6, line7, ...], ...]

    @profile
    def process(self, out_prefix: str, debug_file=None):
        """
        Orchestrates the processing of the GAF file, either in a single thread or in parallel.
        """
        t0 = time.time()
        
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(f"Processing GAF file in parallel with {self.threads} threads...", file=stderr)
        
        # Set global variable before forking to leverage copy-on-write (avoids pickling)
        global shared_align_anchor
        shared_align_anchor = self.align_anchor
        
        # Divide the GAF file into chunks
        gaf_chunks = self._chunk_gaf_file(self.gaf_path, self.threads)
        
        # Initialize the worker processes and run the process_gaf_chunk function on each chunk
        with multiprocessing.Pool(processes=self.threads, initializer=init_worker) as pool:
            results = pool.map(process_gaf_chunk, gaf_chunks)
        
        if settings.DEBUG:
            print("Merging results from worker processes...", file=stderr)
        # Prepare reads_processed TSV: remove old file once before appending
        reads_processed_path = f"{out_prefix}.reads_processed.tsv" if settings.OUTPUT_LOGGING_FILES else None
        if reads_processed_path and os.path.exists(reads_processed_path):
            os.remove(reads_processed_path)
        for result_dict in results:
            self.align_anchor.merge_results(result_dict, reads_processed_path)
        
        total_time_for_gaf_processing = time.time() - t0
        
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(
                f"GAF processing finished in {total_time_for_gaf_processing:.2f}s with {self.threads} threads",
                file=stderr,
            )

        # Run the dump_valid_anchors method which runs the unreliable snarl filtering and the anchor extensions
        
        ########################################
        ## SPECIAL CASE: MIN_ANCHOR_LENGTH = 0
        ########################################
        if settings.MIN_ANCHOR_LENGTH == 0:
            kwargs = {
            "extended_out_file_path": f"{out_prefix}.extended.jsonl",
            "reliable_snarls_out_file_path": f"{out_prefix}.reliable_snarls.tsv",
            "pre_reliable_sizes_out_file_path": f"{out_prefix}.subgraph.sizes.pre_reliable.tsv",
            "path_matched_sizes_out_file_path": f"{out_prefix}.subgraph.sizes.path_matched.tsv",
            }
            if settings.OUTPUT_LOGGING_FILES:
                kwargs.update({
                    "snarl_variant_type_out_file_path": f"{out_prefix}.snarl_variant_type.jsonl",
                    "snarl_compatibility_out_file_path": f"{out_prefix}.snarl_compatibility.jsonl",
                    "snarl_common_reads_out_file_path": f"{out_prefix}.snarl_2_snarl_common_reads.jsonl",
                    "snarl_read_partitions_out_file_path": f"{out_prefix}.snarl_2_snarl_read_partitions.jsonl",
                    "snarl_coverage_out_file_path": f"{out_prefix}.snarl_coverage.jsonl",
                    "snarl_allelic_coverage_out_file_path": f"{out_prefix}.snarl_allelic_coverage.jsonl",
                    "snarl_coverage_extended_out_file_path": f"{out_prefix}.snarl_coverage_extended.jsonl",
                    "snarl_allelic_coverage_extended_out_file_path": f"{out_prefix}.snarl_allelic_coverage_extended.jsonl",
                    "binomial_pairs_out_file_path": f"{out_prefix}.binomial_pairs.tsv"
                })
                self.align_anchor.dump_valid_anchors_0bp(**kwargs)
                self.align_anchor.dump_snarls_and_anchors_in_reads_dict(f"{out_prefix}.snarls_and_anchors_in_reads.jsonl")
            else:
                self.align_anchor.dump_valid_anchors_0bp(**kwargs)
            
            # Always emit the extended subgraph size TSV, independent of OUTPUT_LOGGING_FILES.
            # This is a lightweight summary artifact that downstream steps may rely on.
            out_file = f"{out_prefix}.subgraph.sizes.extended.tsv"
            out_dir = os.path.dirname(out_file)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            self.dump_dict_size_extended(out_file)

        ########################################
        ## NORMAL CASE: MIN_ANCHOR_LENGTH > 0
        ########################################
        if settings.MIN_ANCHOR_LENGTH > 0:
            kwargs = {
                "extended_out_file_path": f"{out_prefix}.extended.jsonl",
                "reliable_snarls_out_file_path": f"{out_prefix}.reliable_snarls.tsv",
                "pre_reliable_sizes_out_file_path": f"{out_prefix}.subgraph.sizes.pre_reliable.tsv",
            }

            if settings.OUTPUT_LOGGING_FILES:
                kwargs.update({
                    "path_matched_sizes_out_file_path": f"{out_prefix}.subgraph.sizes.path_matched.tsv",
                    "seq_matched_sizes_out_file_path": f"{out_prefix}.subgraph.sizes.seq_matched.tsv",
                    "anchor_read_tracking_file_path": f"{out_prefix}.read_drop_tracking.jsonl",
                    "independent_anchor_read_tracking_file_path": f"{out_prefix}.independent_ext_tracking.jsonl",
                    "snarl_variant_type_out_file_path": f"{out_prefix}.snarl_variant_type.jsonl",
                    "snarl_compatibility_out_file_path": f"{out_prefix}.snarl_compatibility.jsonl",
                    "snarl_common_reads_out_file_path": f"{out_prefix}.snarl_2_snarl_common_reads.jsonl",
                    "snarl_read_partitions_out_file_path": f"{out_prefix}.snarl_2_snarl_read_partitions.jsonl",
                    "snarl_coverage_out_file_path": f"{out_prefix}.snarl_coverage.jsonl",
                    "snarl_allelic_coverage_out_file_path": f"{out_prefix}.snarl_allelic_coverage.jsonl",
                    "snarl_coverage_extended_out_file_path": f"{out_prefix}.snarl_coverage_extended.jsonl",
                    "snarl_allelic_coverage_extended_out_file_path": f"{out_prefix}.snarl_allelic_coverage_extended.jsonl",
                    "binomial_pairs_out_file_path": f"{out_prefix}.binomial_pairs.tsv",
                    "adjacent_snarl_pairs_out_file_path": f"{out_prefix}.adjacent_snarl_pairs.jsonl"
                })
                self.align_anchor.dump_valid_anchors(**kwargs)
                self.align_anchor.dump_snarls_and_anchors_in_reads_dict(f"{out_prefix}.snarls_and_anchors_in_reads.jsonl")
            
            else:
                # Just dump the extended valid anchors JSON
                self.align_anchor.dump_valid_anchors(**kwargs)

            # Always emit the extended subgraph size TSV, independent of OUTPUT_LOGGING_FILES.
            # This is a lightweight summary artifact that downstream steps may rely on.
            out_file = f"{out_prefix}.subgraph.sizes.extended.tsv"
            out_dir = os.path.dirname(out_file)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            self.dump_dict_size_extended(out_file)


    def dump_dictionary_with_counts(self, out_file: str):
        """
        It dumps the positioned anchor dictionary by json
        """
        self.align_anchor.dump_dictionary_with_reads_counts(out_file)

    def dump_dict_size_extended(self, out_file: str):
        """
        It dumps the anchors by json
        """
        self.align_anchor.print_extended_anchor_info(out_file) 

    def dump_bandage_csv_extended(self, out_file: str):
        """
        It dumps CSV with node and colour of all anchor nodes
        """
        self.align_anchor.print_sentinels_for_bandage(out_file) 
