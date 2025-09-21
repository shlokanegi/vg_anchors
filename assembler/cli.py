import click
import time
import os.path
from datetime import datetime
import sys
import re

import assembler.constants as constants
from assembler.handler import Orchestrator
from assembler.builder import AnchorDictionary


@click.group()
def cli():
    """vg-anchor is a tool for finding anchors in a variation graph."""
    pass


@cli.command("build", help="Build the anchor dictionary.")
@click.option("--graph",required=True,type=click.Path(exists=True),help="Input packedgraph file (.vg)")
@click.option("--index",required=True,type=click.Path(exists=True),help="Input distance index file (.dist)")
@click.option("--output-prefix", required=True, type=click.Path(), help="Output prefix for the anchor dictionary")
def build(graph, index, output_prefix):
    """Build an anchor dictionary from graph and index files."""
    from assembler.builder import AnchorDictionary

    output_dictionary = output_prefix + ".pkl"
    bandage_csv = output_prefix + ".bandage.csv"
    sizes_csv = output_prefix + ".sizes.tsv"
    
    t0 = time.time()
    dictionary_builder = AnchorDictionary()
    dictionary_builder.build(graph, index)
    dictionary_builder.fill_anchor_dictionary(extend=False)
    print(
        f"Anchors dictionary from {len(dictionary_builder.leaf_snarls)} snarls, containing {len(dictionary_builder.sentinel_to_anchor)} sentinels built in {time.time()-t0:.2f}",
        flush=True,
        file=sys.stderr,
    )
    dictionary_builder.add_positions_to_anchors()
    dictionary_builder.dump_dictionary(output_dictionary)
    dictionary_builder.print_anchor_boundaries_dict(output_prefix)

    if bandage_csv:
        dictionary_builder.print_sentinels_for_bandage(bandage_csv)

    if sizes_csv:
        dictionary_builder.print_dict_sizes(sizes_csv)


@cli.command("get-anchors", help="""Get anchors from a GAF file, given a graph, index, and anchors.""")
@click.option("--dictionary",required=True,type=click.Path(exists=True),help="Input anchor dictionary file")
@click.option("--graph", required=True, type=click.Path(exists=True), help="Input graph file")
@click.option("--alignment",required=True,type=click.Path(exists=True),help="Input alignment file")
@click.option("--fasta",required=True,type=click.Path(exists=True),help="Input fasta file")
@click.option("--output", required=True, type=click.Path(), help="Output basename. Used by anchors (jsonl) and pkl count (.count.pkl)")
@click.option("--threads",default=1,show_default=True,type=click.Path(),help="Number of threads to use for parallel processing.")
def get_anchors(dictionary, graph, alignment, fasta, output, threads):
    """Process alignment and get anchors."""
    from assembler.handler import Orchestrator

    anchors_dir = os.path.dirname(output)
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_path = os.path.join(anchors_dir, "params_run.log")
    log_content = f"""
    VG_ANCHOR PARAMETERS LOG
    Timestamp: {timestamp}
    ==================================================
    
    MIN_ANCHOR_LENGTH = {constants.MIN_ANCHOR_LENGTH}
    EXPECTED_MAP_Q = {constants.EXPECTED_MAP_Q}
    MIN_ANCHOR_READS = {constants.MIN_ANCHOR_READS}
    HET_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING = {constants.HET_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING}
    HOMO_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING = {constants.HOMO_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING}
    MIN_READS_REQUIRED_FOR_MERGING_R0 = {constants.MIN_READS_REQUIRED_FOR_MERGING_R0}
    MIN_READS_REQUIRED_FOR_MERGING_R1 = {constants.MIN_READS_REQUIRED_FOR_MERGING_R1}
    FRACTION_READS_FOR_SNARL_BOUNDARY_EXTENTION = {constants.FRACTION_READS_FOR_SNARL_BOUNDARY_EXTENTION}
    MIN_READS_REQUIRED_FOR_BOUNDARY_EXTENSION = {constants.MIN_READS_REQUIRED_FOR_BOUNDARY_EXTENSION}
    DROP_FRACTION = {constants.DROP_FRACTION}
    MIN_ANCHOR_READCOV = {constants.MIN_ANCHOR_READCOV}

    # PHASING CONSISTENCY CHECK ANCHORS/SNARLS CONSTANTS
    MIN_SNARL_LINKAGE_THRESHOLD = {constants.MIN_SNARL_LINKAGE_THRESHOLD}
    RELIABLE_SNARL_FRACTION_THRESHOLD = {constants.RELIABLE_SNARL_FRACTION_THRESHOLD}
    ADD_BACK_HOMO_SNARLS = {constants.ADD_BACK_HOMO_SNARLS}
    ERROR_TOLERANCE_IN_COMPATIBILITY_CHECK = {constants.ERROR_TOLERANCE_IN_COMPATIBILITY_CHECK}
    ENABLE_UNEQUAL_SET_COMPATIBILITY = {constants.ENABLE_UNEQUAL_SET_COMPATIBILITY}
    MIN_READS_FOR_PARTITION_COMPATIBILITY = {constants.MIN_READS_FOR_PARTITION_COMPATIBILITY}
    """

    with open(log_path, "w") as log_file:
        log_file.write(log_content.strip())

    orchestrator = Orchestrator(dictionary, graph, alignment, fasta, threads)
    orchestrator.process(out_prefix=f"{output}")
    orchestrator.dump_dict_size_extended(f"{output}.subgraph.sizes.extended.tsv")

@cli.command("benchmark-snarl-finding", help="""Benchmark the reliable snarl finding step with multiple thread counts.""")
@click.option("--dictionary",required=True,type=click.Path(exists=True),help="Input anchor dictionary file")
@click.option("--graph", required=True, type=click.Path(exists=True), help="Input graph file")
@click.option("--alignment",required=True,type=click.Path(exists=True),help="Input alignment file")
@click.option("--fasta",required=True,type=click.Path(exists=True),help="Input fasta file")
@click.option("--output", required=True, type=click.Path(), help="Output basename. Used by anchors (jsonl) and pkl count (.count.pkl)")
@click.option("--threads",default=1,show_default=True,type=click.Path(),help="Maximum number of threads to benchmark.")
def benchmark_snarl_finding(dictionary, graph, alignment, fasta, output, threads):
    """Benchmark the reliable snarl finding step."""
    anchors_dir = os.path.dirname(output)
    max_threads = int(threads)
    runtime_logs = []
    with open(os.path.join(anchors_dir, "benchmark_snarl_finding.log"), "w") as log_file:
        print("threads\ttime_for_reliable_snarls_finding", file=log_file)
        log_file.flush()

        for t in range(1, max_threads + 1):
            orchestrator = Orchestrator(dictionary, graph, alignment, fasta, t)
            orchestrator.process(out_prefix=f"{output}")
            log = orchestrator.align_anchor.runtime_logs
            runtime_logs.append(log)
            print(f"{log['threads']}\t{log['time_for_reliable_snarls_finding']:.4f}", file=log_file)
            log_file.flush()
    

if __name__ == "__main__":
    cli()
