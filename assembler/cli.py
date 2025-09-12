import click
import time
import os.path
from datetime import datetime
import sys
import re

from assembler import config


@click.group()
@click.option(
    "--config",
    "config_file",
    type=click.Path(dir_okay=False),
    help="Path to a custom config.ini file. Overrides the default.",
    default=None,
)
def cli(config_file):
    """Anchor processing tool for the assembler package."""
    config.load_config(config_file)


@cli.command()
@click.option(
    "--graph",
    required=True,
    type=click.Path(exists=True),
    help="Input packedgraph file (.vg)",
)
@click.option(
    "--index",
    required=True,
    type=click.Path(exists=True),
    help="Input distance index file (.dist)",
)
@click.option(
    "--output-prefix",
    required=True,
    type=click.Path(),
    help="Output prefix for the anchor dictionary",
)

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
    
    click.echo(f"Anchor dictionary built and saved to {output_dictionary}")


@cli.command()
@click.option(
    "--dictionary",
    required=True,
    type=click.Path(exists=True),
    help="Input anchor dictionary file",
)
@click.option(
    "--graph", required=True, type=click.Path(exists=True), help="Input graph file"
)
@click.option(
    "--alignment",
    required=True,
    type=click.Path(exists=True),
    help="Input alignment file",
)
@click.option(
    "--fasta",
    required=True,
    type=click.Path(exists=True),
    help="Input fasta file"
)
@click.option(
    "--output", required=True, type=click.Path(), help="Output basename. Used by anchors (jsonl) and pkl count (.count.pkl)"
)
def get_anchors(dictionary, graph, alignment, fasta, output):
    """Process alignment and get anchors."""
    from assembler.handler import Orchestrator

    anchors_dir = os.path.dirname(output)
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_path = os.path.join(anchors_dir, "params_run.log")
    log_content = f"""
    VG_ANCHOR PARAMETERS LOG
    Timestamp: {timestamp}
    ==================================================
    
    MIN_ANCHOR_LENGTH = {config.settings.getint('MIN_ANCHOR_LENGTH')}
    EXPECTED_MAP_Q = {config.settings.getint('EXPECTED_MAP_Q')}
    MIN_ANCHOR_READS = {config.settings.getint('MIN_ANCHOR_READS')}
    HET_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING = {config.settings.getfloat('HET_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING')}
    HOMO_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING = {config.settings.getfloat('HOMO_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING')}
    MIN_READS_REQUIRED_FOR_MERGING_R0 = {config.settings.getint('MIN_READS_REQUIRED_FOR_MERGING_R0')}
    MIN_READS_REQUIRED_FOR_MERGING_R1 = {config.settings.getint('MIN_READS_REQUIRED_FOR_MERGING_R1')}
    FRACTION_READS_FOR_SNARL_BOUNDARY_EXTENTION = {config.settings.getfloat('FRACTION_READS_FOR_SNARL_BOUNDARY_EXTENTION')}
    MIN_READS_REQUIRED_FOR_BOUNDARY_EXTENSION = {config.settings.getint('MIN_READS_REQUIRED_FOR_BOUNDARY_EXTENSION')}
    DROP_FRACTION = {config.settings.getfloat('DROP_FRACTION')}
    MIN_ANCHOR_READCOV = {config.settings.getint('MIN_ANCHOR_READCOV')}

    # PHASING CONSISTENCY CHECK ANCHORS/SNARLS CONSTANTS
    MIN_SNARL_LINKAGE_THRESHOLD = {config.settings.getint('MIN_SNARL_LINKAGE_THRESHOLD')}
    RELIABLE_SNARL_FRACTION_THRESHOLD = {config.settings.getfloat('RELIABLE_SNARL_FRACTION_THRESHOLD')}
    ADD_BACK_HOMO_SNARLS = {config.settings.getboolean('ADD_BACK_HOMO_SNARLS')}
    ERROR_TOLERANCE_IN_COMPATIBILITY_CHECK = {config.settings.getint('ERROR_TOLERANCE_IN_COMPATIBILITY_CHECK')}
    ENABLE_UNEQUAL_SET_COMPATIBILITY = {config.settings.getboolean('ENABLE_UNEQUAL_SET_COMPATIBILITY')}
    MIN_READS_FOR_PARTITION_COMPATIBILITY = {config.settings.getint('MIN_READS_FOR_PARTITION_COMPATIBILITY')}
    """

    with open(log_path, "w") as log_file:
        log_file.write(log_content.strip())

    t1 = time.time()
    orchestrator = Orchestrator(dictionary, graph, alignment, fasta)
    orchestrator.process(f"{output}")
    print(
        f"GAF alignment processed in {time.time()-t1:.2f}", flush=True, file=sys.stderr
    )

    orchestrator.dump_anchors(f"{output}.jsonl", f"{output}.extended.jsonl", f"{output}.anchor_reads_tracker.jsonl", f"{output}.independent_extension.jsonl", f"{output}.extended.pruned.jsonl", f"{output}.reliable_snarls.tsv", f"{output}.snarl_variant_type.jsonl", f"{output}.snarl_compatibility.jsonl", f"{output}.snarl_2_snarl_common_reads.jsonl", f"{output}.snarl_2_snarl_read_partitions.jsonl", f"{output}.snarl_coverage.jsonl", f"{output}.snarl_allelic_coverage.jsonl", f"{output}.snarl_coverage_extended.jsonl", f"{output}.snarl_allelic_coverage_extended.jsonl")
    orchestrator.dump_dict_size_extended(f"{output}.subgraph.sizes.extended.tsv")

if __name__ == "__main__":
    cli()
