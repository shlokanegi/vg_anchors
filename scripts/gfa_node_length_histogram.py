from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.axis as ax
import matplotlib.gridspec as gridspec
import seaborn as sns
import sys
import numpy as np
from scipy import stats

def get_len_dictionary_from_gfa(filename) -> dict :
    node_length_dictionary = defaultdict(int)
    with open(filename, "r") as f:
        for line in f:
            line_stripped = line.strip()
            if line_stripped.startswith('S'):
                _, node_name, seq = line_stripped.split('\t')
                seq_len = len(seq)
                if seq_len > 1:
                    node_length_dictionary[seq_len] += 1
                if seq_len > 10000:
                    print(node_name)
    return node_length_dictionary

def node_length_distribution(filename1: str, filename2: str, outfile: str) -> None:

    node_lengths1 = get_len_dictionary_from_gfa(filename1)
    node_lengths2 = get_len_dictionary_from_gfa(filename2)

    title1 = "Original Dataset Distribution"
    title2 = "VG Modified Dataset Distribution"

    plt.figure(figsize=(12, 12))
    gs = gridspec.GridSpec(2, 2)
    ax1 = plt.subplot(gs[0, :])
    ax2 = plt.subplot(gs[1, :-1])
    ax3 = plt.subplot(gs[1, 1:])

    get_scatter_plot(ax1, node_lengths1, node_lengths2)
    get_violin_plot(ax2, node_lengths1, title1,log_scale=(False,True))
    get_violin_plot(ax3, node_lengths2, title2)

    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

def get_scatter_plot(axis: ax, dataset1: dict, dataset2: dict) -> None:
    axis.scatter(list(dataset1.keys()), list(dataset1.values()), 
                      label='Original', alpha=0.7, color='blue')
    axis.scatter(list(dataset2.keys()), list(dataset2.values()), 
                      label='Adjusted by vg', alpha=0.7, color='red')
    axis.set_yscale('log')
    axis.set_title('Node length dot plot')
    axis.set_xlabel('Node Length')
    axis.set_ylabel('Frequency')
    axis.legend()
    
def get_violin_plot(axis: ax, dataset: dict, title: str, log_scale=None) -> None:

    df_data = []
    for val, count in dataset.items():
        df_data.extend([val] * count)
    
    # Violin plot with more details
    sns.violinplot(y=df_data, ax=axis, inner='box', cut=0, log_scale=log_scale)
    axis.set_title(title)
    axis.set_xlabel('Node Length')
    axis.set_ylabel('Value')
    
    # Add some statistical annotations
    axis.text(0.05, 0.95, f'Mean: {np.mean(df_data):.2f}\nMedian: {np.median(df_data):.2f}', 
             transform=axis.transAxes, verticalalignment='center',horizontalalignment='right',)


if __name__ == "__main__":
    node_length_distribution(sys.argv[1], sys.argv[2], sys.argv[3])
