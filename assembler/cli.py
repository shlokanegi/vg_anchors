import click
import time
import os.path
import sys

from assembler.constants import MIN_ANCHOR_LENGTH
from assembler.handler import Orchestrator
from assembler.builder import AnchorDictionary
import assembler.qc
import assembler.helpers


@click.group()
def cli():
    """Anchor processing tool for the assembler package."""
    pass


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
# @click.option("--anchors-json", type=click.Path(), help="Output file for the anchors in the dictionary (.json)")
# @click.option("--bandage-csv", type=click.Path(), help="Output CSV file for Bandage")
# @click.option("--sizes-csv", type=click.Path(), help="Output CSV file for anchor sizes")
# @click.option(
#     "--positioned-dict", type=click.Path(), help="Output file for positioned dictionary"
# )
def build(graph, index, output_prefix):
    output_dictionary = output_prefix + ".pkl"
    bandage_csv = output_prefix + ".bandage.csv"
    sizes_csv = output_prefix + ".sizes.tsv"
    # paths_file = output_prefix + ".used_pathnames.txt"
    # positioned_dict = output_prefix + ".positioned.json"

    """Build an anchor dictionary from graph and index files."""
    t0 = time.time()
    dictionary_builder = AnchorDictionary()
    dictionary_builder.build(graph, index)
    dictionary_builder.fill_anchor_dictionary()
    print(
        f"Anchors dictionary from {len(dictionary_builder.leaf_snarls)} snarls, containing {len(dictionary_builder.sentinel_to_anchor)} sentinels built in {time.time()-t0:.2f}",
        flush=True,
        file=sys.stderr,
    )
    dictionary_builder.add_positions_to_anchors()
    dictionary_builder.dump_dictionary(output_dictionary)

    if bandage_csv:
        dictionary_builder.print_sentinels_for_bandage(bandage_csv)

    if sizes_csv:
        dictionary_builder.print_dict_sizes(sizes_csv)
    
    # if paths_file:
    #     dictionary_builder.print_paths_used(paths_file)

    # if positioned_dict:
    #     dictionary_builder.generate_positioned_dictionary("", positioned_dict)

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
    "--output", required=True, type=click.Path(), help="Output file for anchors"
)
def get_anchors(dictionary, graph, alignment, output):
    # positioned_dict = dictionary.rstrip("pkl") + "positioned.json"
    """Process alignment and get anchors."""
    t1 = time.time()
    orchestrator = Orchestrator(dictionary, graph, alignment)
    orchestrator.process(f"{output}.anchors_info.csv")
    print(
        f"GAF alignment processed in {time.time()-t1:.2f}", flush=True, file=sys.stderr
    )

    orchestrator.dump_anchors(output)
    orchestrator.dump_dictionary_with_counts(dictionary.rstrip("pkl") + "count.pkl")
    click.echo(f"Anchors have minimun length of {MIN_ANCHOR_LENGTH}")
    click.echo(f"Anchors processed and saved to {output}")

@cli.command()
@click.option(
    "--anchors",
    required=True,
    type=click.Path(exists=True),
    help="Input anchors obtained using get_anchors",
)
@click.argument(
    "fastq", 
    required=True,
    nargs=-1,  # Allow multiple fastq files as arguments
    type=click.Path(exists=True)
)
def verify_output(anchors, fastq):

    anchors_name = anchors.split('/')[-1].split('.')[0]
    
    fastq_stripped = fastq[0].rstrip(".fastq") if fastq[0].endswith(".fastq") else fastq[0].rstrip(".fastq.gz")
    fastq_name = fastq_stripped.split('/')[-1]
    fastq_path = fastq_stripped.rstrip(fastq_name)
    out_fastq = fastq_path + f"{anchors_name}.selected.fastq"
    assembler.qc.verify_anchors_validity(anchors, fastq, out_fastq)

@cli.command()
# @click.option(
#     "--anchors-dict",
#     required=True,
#     type=click.Path(exists=True),
#     help="Input anchors computed",
# )
@click.option(
    "--anchors-count",
    required=True,
    type=click.Path(exists=True),
    help="Input anchors count ",
)
@click.option(
    "--plot-title",
    required=True,
    help="Title of the plot ",
)
@click.option(
    "--out-png", required=True, help="prefix of the png files in output"
)
def plot_stats( anchors_count, out_png, plot_title):

    assembler.helpers.plot_count_histogram(anchors_count, out_png + "count.png")

    assembler.helpers.plot_anchor_count_genome_distribution(
        anchors_count, out_png + "position_count.png", plot_title,
    )
    assembler.helpers.plot_heteroxigosity_on_genome(anchors_count, out_png + "het.png", plot_title)
    



if __name__ == "__main__":
    cli()
