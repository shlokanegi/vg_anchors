import click
import time
import os
from datetime import datetime
import sys
import re

from assembler.config import settings, load_config
from assembler.handler import Orchestrator
from assembler.builder import AnchorDictionary

# Debugpy setup - only enabled if VG_ANCHORS_DEBUG environment variable is set
# if os.environ.get("VG_ANCHORS_DEBUG", "").lower() in ("1", "true", "yes"):
#     import debugpy
#     port = int(os.environ.get("VG_ANCHORS_DEBUG_PORT", "5681"))
#     debugpy.listen(('127.0.0.1', port))
#     print(f"Waiting for debugger to attach on port {port}...", flush=True, file=sys.stderr)
#     debugpy.wait_for_client()
#     print(f"Debugger attached! Resuming script.", flush=True, file=sys.stderr)

@click.group()
@click.option(
    "--config",
    "config_file",
    type=click.Path(dir_okay=False),
    help="Path to a custom config.ini file. Overrides the default.",
    default=None,
)
def cli(config_file):
    """vg-anchor is a tool for finding anchors in a variation graph."""
    load_config(config_file)


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
@click.option("--shasta2", is_flag=True, help="Enable Shasta2-specific anchor generation.")
def get_anchors(dictionary, graph, alignment, fasta, output, threads, shasta2):
    """Process alignment and get anchors."""
    from assembler.handler import Orchestrator

    read_id_map = None
    if shasta2:
        read_id_map = {}
        with open(fasta, "r") as f:
            read_id = 0
            for line in f:
                if line.startswith(">"):
                    read_name = line.strip().split()[0][1:]
                    read_id_map[read_name] = read_id
                    read_id += 1

    anchors_dir = os.path.dirname(output)
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_path = os.path.join(anchors_dir, "params_run.log")

    # Dynamically generate the constants log
    constants_log_lines = []
    if settings.raw_config:
        for section in settings.raw_config.sections():
            constants_log_lines.append(f"[{section}]")
            for key, value in settings.raw_config.items(section):
                constants_log_lines.append(f"{key.upper()} = {value}")
    
    log_content = f"""
    VG_ANCHOR PARAMETERS LOG
    Timestamp: {timestamp}
    ==================================================
    
    """ + "\n    ".join(constants_log_lines) + "\n"

    with open(log_path, "w") as log_file:
        log_file.write(log_content)

    orchestrator = Orchestrator(
        dictionary_path=dictionary,
        graph_path=graph,
        gaf_path=alignment,
        fasta_path=fasta,
        threads=threads,
        read_id_map=read_id_map
    )
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
    

@cli.command("chunker-build", help="Build the chunker dictionary.")
@click.option("--graph",required=True,type=click.Path(exists=True),help="Input packedgraph file (.vg)")
@click.option("--index",required=True,type=click.Path(exists=True),help="Input distance index file (.dist)")
@click.option("--output-prefix", required=True, type=click.Path(), help="Output prefix for the chunker dictionary")
def chunker_build(graph, index, output_prefix):
    """Build an anchor dictionary from graph and index files."""
    from assembler.chunker_builder import ChunkerAnchorDictionary

    hap_counts_tsv = output_prefix + ".hap_counts.tsv"
    
    t0 = time.time()
    dictionary_builder = ChunkerAnchorDictionary()
    dictionary_builder.build(graph, index)
    
    # Process leaf snarls and generate boundaries (required before step counting)
    print("Processing leaf snarls...", flush=True, file=sys.stderr)
    dictionary_builder.process_snarls()
    print(f"Found {len(dictionary_builder.leaf_snarls)} leaf snarls", flush=True, file=sys.stderr)
    
    print("Generating snarl boundaries...", flush=True, file=sys.stderr)
    dictionary_builder.generate_anchors_boundaries()
    print(f"Processed {len(dictionary_builder.snarl_boundaries[0])} snarl boundaries", flush=True, file=sys.stderr)
    
    dictionary_builder.get_step_counts_from_sentinel_nodes_of_snarls()
    print(
        f"Step counts from sentinel nodes of snarls built in {time.time()-t0:.2f}",
        flush=True,
        file=sys.stderr,
    )
    
    dictionary_builder.print_hap_counts(hap_counts_tsv)

if __name__ == "__main__":
    cli()
