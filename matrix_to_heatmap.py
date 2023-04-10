#! /usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt

# parses command-line arguments with coverage tables and read counts
def parse_args():
    parser = argparse.ArgumentParser(description="Generates heatmap from count matrix")
    parser.add_argument("counts", help="tsv in the form feature \t count_1 \t ... count_x")
    parser.add_argument("--out", "-o", default="heatmap", help="prefix for heatmap")
    return parser.parse_args()

def parse_read_counts(count_matrix, output):
    df = pd.read_csv(count_matrix, sep=",|\t", index_col=0, engine="python")
    # Displaying dataframe as an heatmap
    # with diverging colourmap as RdYlBu
    plt.figure(figsize=(10, 60))
    im = plt.imshow(df, cmap ="RdYlBu_r")
    
    # Displaying a color bar to understand
    # which color represents which range of data
    plt.colorbar()
    
    # Assigning labels of x-axis 
    # according to dataframe
    plt.xticks(range(len(df.columns)), df.columns, rotation = 315, ha="left", va="top", rotation_mode="anchor", fontsize=15)
    
    # Assigning labels of y-axis 
    # according to dataframe
    plt.yticks(range(len(df.index)), df.index, fontsize=15)

    bottom, top = plt.ylim()
    plt.ylim(top=top-0.5)  # adjust the top leaving bottom unchanged
    plt.ylim(bottom=bottom+0.5)  # adjust the bottom leaving top unchanged
    #ax.set_ylim(bottom + 0.5, top - 0.5)
    
    #plt.setp(plt.xticks, rotation=45, ha="left", va="top", rotation_mode="anchor", fontsize=10)

    # Displaying the figure
    plt.savefig(output + ".jpg", bbox_inches="tight")
    #plt.show()
    """
    # generates a heatmap based on orgaism similarity
    fig, ax = plt.subplots()

    #ax.set_xticks(range(len(df.columns)), df.columns); ax.set_yticks(range(len(df.index)), df.index)
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels(df.columns)
    ax.set_yticks(range(len(df.index)))
    if len(df.index) < 30:
        ax.set_yticklabels(df.index)
    #ax.set_title("{0}-mer Distance Estimation".format(k))

    # creates heatmap from axes and matrix
    ax.xaxis.tick_top()
    heatmap = ax.imshow(df, cmap ="RdYlBu_r")
    fig.colorbar(heatmap)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="left", 
        va="top", rotation_mode="anchor", fontsize=10)

    plt.setp(ax.get_yticklabels(), fontsize=10)

    #fig.set_size_inches(10, 8)
    plt.savefig(output + ".jpg", bbox_inches="tight")
    print("Heatmap saved to {0}.jpg".format(output))
    """



def main():
    args = parse_args()
    parse_read_counts(args.counts, args.out)
    
main()