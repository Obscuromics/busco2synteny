#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: busco2synteny.py -a <STR> -b <STR> -x <STR> -y <STR> [-l <STR> -m <INT> -n <INT> -c <INT> -t <FLT> -w <FLT> -h]

  [Options]
    -a, --genomefileA <STR>                     Genomefile for taxon A
    -b, --genomefileB <STR>                     Genomefile for taxon B
    -x, --buscoA <STR>                          BUSCO 'full_table.tsv' for taxon A
    -y, --buscoB <STR>                          BUSCO 'full_table.tsv' for taxon B
    -l, --labels <STR>                          Whether to plot labels, choose from True or False [default: True]
    -m, --gapA <INT>                            Gap between genome A chromosomes [default: 10_000_000]
    -n, --gapB <INT>                            Gap between genome B chromosomes [default: 10_000_000]
    -c, --chromosome_width <INT>                Chromosome width [default: 6]
    -t, --alpha <FLT>                           Alpha of alignments [default: 0.1]
    -w, --linewidth <FLT>                           Linewidth of alignments [default: 0.1]
    -h, --help                                  Show this message

"""

import sys
from docopt import docopt
import collections
from scipy.interpolate import make_interp_spline, BSpline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from operator import itemgetter

# mamba install docopt scipy numpy pandas matplotlib pyarrow

# Modified version of a script originally created by Alex Mackintosh
# https://github.com/A-J-F-Mackintosh/Mackintosh_et_al_2022_Binodaphne/blob/main/minimap2synteny.py


def generate_genomefile_dict(genomefile, offset, colour):
    genomefile_dict = {}
    cumulative_genome = offset
    with open(genomefile, "r") as fin:
        # for each chromosome, record cumulative coordinates, orientation, and label
        for i, line in enumerate(fin):
            line = line.rstrip()
            chromosome, chromosome_length, orientation, label = line.split("\t")
            chromosome_length = int(chromosome_length)
            genomefile_dict[chromosome] = [
                cumulative_genome,
                cumulative_genome + chromosome_length,
                orientation,
                label,
            ]
            cumulative_genome += chromosome_length
            cumulative_genome += offset
            if colour:
                genomefile_dict[chromosome].append(i)
    return genomefile_dict


# plots chromosomes as lines and adds labels if arg is True
def plot_chromosomes(genomefile_dict, y_coord, labels, chromosome_width):
    for chromosome in genomefile_dict:
        plt.plot(
            [genomefile_dict[chromosome][0], genomefile_dict[chromosome][1]],
            [y_coord, y_coord],
            color="slategrey",
            alpha=1,
            linewidth=chromosome_width,
        )
        middle_of_chromosome = genomefile_dict[chromosome][0] + (
            (genomefile_dict[chromosome][1] - genomefile_dict[chromosome][0]) / 2
        )
        plt.plot(
            [genomefile_dict[chromosome][0], genomefile_dict[chromosome][1]],
            [y_coord * 1.02, y_coord * 1.02],
            color="slategrey",
            alpha=1,
            linewidth=chromosome_width,
        )
        middle_of_chromosome = genomefile_dict[chromosome][0] + (
            (genomefile_dict[chromosome][1] - genomefile_dict[chromosome][0]) / 2
        )
        if labels == "True":
            plt.text(
                middle_of_chromosome,
                y_coord * 1.2,
                genomefile_dict[chromosome][3],
                ha="center",
                va="center",
                wrap=True,
                fontsize=20,
            )


def generate_alignment_dicts(
    liftover, genomefile_A_dict, genomefile_B_dict, alignment_coords
):
    with open(liftover, "r") as fin:
        # for each alignment, record coordinates in either genome
        for line in fin:
            alignment = []
            line = line.rstrip()
            seqA, startA, endA, seqB, startB, endB, seqanc = line.split("\t")
            midpointA = int(startA) + ((int(endA) - int(startA)) / 2)
            midpointB = int(startB) + ((int(endB) - int(startB)) / 2)
            average_width = ((int(endA) - int(startA)) + (int(endB) - int(startB))) / 2

            if seqA in genomefile_A_dict.keys():
                # flip coords if orientation is -
                if genomefile_A_dict[seqA][2] == "-":
                    midpointA = genomefile_A_dict[seqA][1] - midpointA
                if genomefile_A_dict[seqA][2] == "+":
                    midpointA += genomefile_A_dict[seqA][0]
                alignment.append(midpointA)

            if seqB in genomefile_B_dict.keys():
                # flip coords if orientation is -
                if genomefile_B_dict[seqB][2] == "-":
                    midpointB = genomefile_B_dict[seqB][1] - midpointB
                if genomefile_B_dict[seqB][2] == "+":
                    midpointB += genomefile_B_dict[seqB][0]
                alignment.append(midpointB)
                if int(seqanc) == 0:
                    line_colour = "lightgrey"
                else:
                    line_colour = cm.tab20((int(seqanc) - 1) / 7)

            # only interested in alignments on sequences in both genomefiles
            if len(alignment) == 2:
                alignment.append(average_width)
                alignment.append(line_colour)
                alignment_coords.append(alignment)

    return alignment_coords


# code from https://stackoverflow.com/questions/19394505/expand-the-line-with-specified-width-in-data-unit/42972469#42972469
def linewidth_from_data_units(linewidth, axis, reference="x"):
    """
    Convert a linewidth in data units to linewidth in points.

    Parameters
    ----------
    linewidth: float
            Linewidth in data units of the respective reference-axis
    axis: matplotlib axis
            The axis which is used to extract the relevant transformation
            data (data limits and size must not change afterwards)
    reference: string
            The axis that is taken as a reference for the data width.
            Possible values: 'x' and 'y'. Defaults to 'y'.

    Returns
    -------
    linewidth: float
            Linewidth in points
    """
    fig = axis.get_figure()
    if reference == "x":
        length = fig.bbox_inches.width * axis.get_position().width
        value_range = np.diff(axis.get_xlim())
    elif reference == "y":
        length = fig.bbox_inches.height * axis.get_position().height
        value_range = np.diff(axis.get_ylim())
    # Convert length to points
    length *= 72
    # Scale linewidth to value range
    return linewidth * (length / value_range)


def load_busco_results(file):
    busco_colnames = ["busco_id", "status", "seq_code", "start", "stop"]

    df = pd.read_csv(
        file,
        sep="\t",
        usecols=[0, 1, 2, 3, 4],
        header=None,
        names=busco_colnames,
        skiprows=[0, 1, 2],
    )

    df = df[df["status"] == "Complete"]
    df.drop(labels=["status"], axis=1, inplace=True)
    if "|" in df["seq_code"].iloc[0]:
        df["seq"] = (
            df["seq_code"]
            .str.split("|")
            .apply(itemgetter(1))
            .str.split(":")
            .apply(itemgetter(0))
        )
    else:
        df["seq"] = df["seq_code"]
    df.drop(labels=["seq_code"], axis=1, inplace=True)
    return df


def get_labels(seqs):
    labels = []
    for i in range(len(seqs)):
        labels.append(int(i + 1))
    return labels


def label_colours(df):
    seqs = sorted(df["seq_A"].unique())
    labels = get_labels(seqs)
    for seq, label in zip(seqs, labels):
        df.loc[df["seq_A"] == seq, "colour"] = label


def create_liftover_from_busco(file_1, file_2):
    df1 = load_busco_results(file_1)
    df2 = load_busco_results(file_2)

    df3 = df1.merge(df2, how="outer", on="busco_id")
    df3.dropna(inplace=True)
    df3.rename(
        columns={
            "start_x": "start_A",
            "stop_x": "end_A",
            "seq_x": "seq_A",
            "start_y": "start_B",
            "stop_y": "end_B",
            "seq_y": "seq_B",
        },
        inplace=True,
    )

    label_colours(df3)
    df3["start_A"] = df3.start_A.astype(int)
    df3["end_A"] = df3.end_A.astype(int)
    df3["start_B"] = df3.start_B.astype(int)
    df3["end_B"] = df3.end_B.astype(int)
    df3["colour"] = df3.colour.astype(int)

    df3[["seq_A", "start_A", "end_A", "seq_B", "start_B", "end_B", "colour"]].to_csv(
        "liftover.tsv", sep="\t", index=False, header=False
    )

if __name__ == "__main__":

    args = docopt(__doc__)

    # generate liftover file
    create_liftover_from_busco(args["--buscoA"], args["--buscoB"])

    # generate dicts for each genome with cumulative coordinates
    genomefile_A_dict = generate_genomefile_dict(
        args["--genomefileA"], int(args["--gapA"]), colour=False
    )
    genomefile_B_dict = generate_genomefile_dict(
        args["--genomefileB"], int(args["--gapB"]), colour=True
    )

    # set up plot
    fig = plt.figure(figsize=(28, 4), frameon=False)
    ax = fig.add_subplot(111)
    ax.axis("off")
    plt.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis="y", which="both", right=False, left=False, labelleft=False)

    # each alignment has coordinates to be recorded
    alignment_coords = []
    alignment_coords = generate_alignment_dicts(
        "liftover.tsv", genomefile_A_dict, genomefile_B_dict, alignment_coords
    )

    # for each alignment, if the two entries are from the two genomes, plot a curve between them
    for alignment in alignment_coords:
        # code is adapted from https://gist.github.com/andrewgiessel/5684769
        upper_x = alignment[0]
        lower_x = alignment[1]
        average_width = alignment[2]
        alignment_colour = alignment[3]
        Y = np.array([-1, -0.568, -0.32, -0.16, -0.056, 0, 0.056, 0.16, 0.32, 0.568, 1])
        X = np.array(
            [
                lower_x,
                lower_x + (0.1 * (upper_x - lower_x)),
                lower_x + (0.2 * (upper_x - lower_x)),
                lower_x + (0.3 * (upper_x - lower_x)),
                lower_x + (0.4 * (upper_x - lower_x)),
                lower_x + (0.5 * (upper_x - lower_x)),
                lower_x + (0.6 * (upper_x - lower_x)),
                lower_x + (0.7 * (upper_x - lower_x)),
                lower_x + (0.8 * (upper_x - lower_x)),
                lower_x + (0.9 * (upper_x - lower_x)),
                upper_x,
            ]
        )
        # requires sorted arrays, so flip if in the wrong order
        if lower_x > upper_x:
            X = np.flip(X)
            Y = np.flip(Y)
        xnew = np.linspace(X.min(), X.max(), 300)
        spl = make_interp_spline(X, Y, k=3)
        power_smooth = spl(xnew)
        linewidth = linewidth_from_data_units(average_width, ax, reference="x")

        plt.plot(
            xnew,
            power_smooth,
            color=alignment_colour,
            alpha=float(args["--alpha"]),
            linewidth=float(args["--linewidth"])
        )
        

    # plot the chromosomes
    plot_chromosomes(
        genomefile_A_dict, 1, args["--labels"], int(args["--chromosome_width"])
    )
    plot_chromosomes(
        genomefile_B_dict, -1, args["--labels"], int(args["--chromosome_width"])
    )

    plt.savefig("busco2synteny.pdf", format="pdf", bbox_inches="tight")

"""
# example command
busco2synteny.py -a fol_ang.genomefile.tsv -b smi_aqu.genomefile.tsv -x busco_results/fol_ang/run_arthropoda_odb10/full_table.tsv -y busco_results/smi_aqu/run_arthropoda_odb10/full_table.tsv -t 0.25
"""
