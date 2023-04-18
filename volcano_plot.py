#! /usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# parses command-line arguments 
def parse_args():
    parser = argparse.ArgumentParser(description="Creates volcano plot from deseq results")
    parser.add_argument("matrix", help="Deseq results with padj and log2FoldChange columns")
    parser.add_argument("--alpha", "-p", type=float, default=0.05, help="P-value and adjusted p-value significance threshold")
    parser.add_argument("--lfc", "-l", type=float, default=1, help="Log fold change threshold")
    parser.add_argument("--valpha", "-v", type=float, default=0.4, help="Point transparency")
    parser.add_argument("--out", "-o", help="output prefix for all files")
    return parser.parse_args()

def make_plot(counts_matrix, lfc_thr=1, pv_thr=0.05, valpha=0.4, figname="volcano"):
    df = pd.read_csv(counts_matrix, sep=",|\t", index_col=0, engine="python")

    # convert data into volcano plot arrays
    x = df['log2FoldChange'].to_numpy()
    y = np.log10(df['padj'].to_numpy())
    y = -1*y

    # color logic
    colors = []
    for idx, yval in enumerate(y):
        if yval < -1*np.log10(pv_thr):
            colors.append("gray")
            continue
        if x[idx] > lfc_thr:
            colors.append("green")
            continue
        if x[idx] < -1*lfc_thr:
            colors.append("red")
            continue
        colors.append("gray")

    # scatter plot
    fig, ax = plt.subplots()
    ax.scatter(x, y, c=colors, alpha=valpha)
    # label axes
    ax.set_xlabel("log2FoldChange")
    ax.set_ylabel("-log10(P)")

    # get x axis bounds
    low, high = plt.xlim()
    xbound = max(abs(low), abs(high))
    plt.xlim(-xbound, xbound)

    # get y axis bounds
    low, high = plt.ylim()

    # draw volcano plot lines
    ax.hlines(y=-np.log10(pv_thr), xmin=-xbound, xmax=xbound, color="gray", linestyle="dashed")
    ax.vlines(x=-lfc_thr, ymin=low, ymax=high, color="gray", linestyle="dashed")
    ax.vlines(x=lfc_thr, ymin=low, ymax=high, color="gray", linestyle="dashed")

    plt.savefig("{0}.png".format(figname))


def main():
    args = parse_args()
    make_plot(args.matrix, args.lfc, args.alpha, args.valpha, args.out)

main()
