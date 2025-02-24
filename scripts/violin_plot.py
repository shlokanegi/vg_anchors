import sys
import csv
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import pandas as pd
import numpy as np

csv.field_size_limit(sys.maxsize)

MAPQ_COL = 11
SEQDIV_COL = 15
GAF_FIELDS = ["Query name","Querylength","Query start","Query end","Strand","Path","Path length","Path start","Path end","# Residue matches","Alignment_block_length","MapQ","DP Alignment Score","bq","Difference sequence","Approx per-base sequence divergence"]
#dv:f:0.0003
def create_violin_plot(tsv_file, output_file):

    tsv_file_name = tsv_file.strip().split('/')[-1].split(".")[0]

    with open(tsv_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        data = list(reader)

    x_column = []
    y_column = []
    for row in data:
        if len(row) > len(GAF_FIELDS):
            x_column.append(int(row[MAPQ_COL + 2]))
            y_column.append(float(row[SEQDIV_COL + 2][5:]))  # Convert y values to float
        else:
            x_column.append(int(row[MAPQ_COL]))
            y_column.append(float(row[SEQDIV_COL][5:]))  # Convert y values to float

    x_counts = Counter(x_column)
    x_unique = sorted(set(x_column))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 12), sharex=True)

    # Violin Plot
    violin_parts = ax1.violinplot(
        [np.array([y for x, y in zip(x_column, y_column) if x == pos]) for pos in x_unique],
        positions=x_unique,
        showmeans=False,
        showextrema=False,
        showmedians=True
    )
    for pc in violin_parts['bodies']:
        pc.set_facecolor('#D43F3A')
        pc.set_edgecolor('black')
        pc.set_alpha(0.7)

    x_60 = [60] * len([y for x, y in zip(x_column, y_column) if x == 60])
    y_60 = [y for x, y in zip(x_column, y_column) if x == 60]
    ax1.scatter(x_60, y_60, color='black', alpha=0.5, s=10)

    ax1.set_title(f'{GAF_FIELDS[SEQDIV_COL]} Distribution by {GAF_FIELDS[MAPQ_COL]} for {tsv_file_name}')
    ax1.set_ylabel(GAF_FIELDS[SEQDIV_COL])

    # Histogram
    ax2.bar(x_unique, [x_counts[x] for x in x_unique], align='center')
    ax2.set_title(f'Histogram: Count of rows for each {GAF_FIELDS[MAPQ_COL]}')
    ax2.set_xlabel(GAF_FIELDS[MAPQ_COL])
    ax2.set_ylabel('Count')

    # Set x-ticks for both plots
    plt.xticks(x_unique)

    # Rotate x-axis labels for both plots
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # Adjust subplot spacing
    plt.subplots_adjust(hspace=0.3)

    # Ensure x-axis labels are visible
    fig.autofmt_xdate()

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

if __name__ == "__main__":
    # Get user input
    tsv_file = sys.argv[1]
    out_png = sys.argv[2]

    # Create the violin plot
    create_violin_plot(tsv_file, out_png)