#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: busco3synteny.py -a <STR> (-x <STR> | -y <STR>)  [-l <STR> -m <INT> -n <INT> -c <INT> -t <FLT> (-w <FLT> | -s) -h]

  [Options]
    -a, --genomefiles <STR>                     File of genome file paths
    -x, --buscofiles <STR>                      File of busco file paths
    -y, --liftoverfiles <STR>                   File of liftover file paths
    -l, --labels <STR>                          Whether to plot labels, choose from True or False [default: True]
    -m, --gapA <INT>                            Gap between genome A chromosomes [default: 10_000_000]
    -n, --gapB <INT>                            Gap between genome B chromosomes [default: 10_000_000]
    -c, --chromosome_width <INT>                Chromosome width [default: 6]
    -t, --alpha <FLT>                           Alpha of alignments [default: 0.1]
    -w, --linewidth <FLT>                       Linewidth of alignments [default: 1]
    -s, --scaled                                Scale linewidths by aligned region length
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

# Generalising busco2synteny.py for N genomes
# Modified version of a script originally created by Alex Mackintosh
# https://github.com/A-J-F-Mackintosh/Mackintosh_et_al_2022_Binodaphne/blob/main/minimap2synteny.py

# genomefile example
"""
SEQ1    56704976    +    Fol_ang_1
SEQ2    24513281    +    Fol_ang_2
SEQ3    14204106    +    Fol_ang_3
SEQ4    13484354    +    Fol_ang_4
SEQ5    9182351    +    Fol_ang_5
"""

"""
How to do consistent line colouring?
    Instead of generating a liftover per pair, might be worth trying to make an overall liftover file
    Can then label colours by the number of genomes they are in per chromosome
Include all buscos but will have to figure out how to do duplicates
Clear nametags for intermediate genomes?
"""


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
            try:
                seqA, startA, endA, seqB, startB, endB, seqanc = line.split("\t")
            except ValueError:
                next
            midpointA = int(float(startA)) + ((int(endA) - int(startA)) / 2)
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

def plot_pair(
    i,
    genomefile_A_f,
    genomefile_B_f,
    busco_A_f=None,
    busco_B_f=None,
    liftover_f="liftover.tsv",
):
    # generate liftover file if using BUSCOs
    if busco_A_f:
        create_liftover_from_busco(busco_A_f, busco_B_f)

    # generate dicts for each genome with cumulative coordinates
    genomefile_A_dict = generate_genomefile_dict(
        genomefile_A_f, int(args["--gapA"]), colour=False
    )
    genomefile_B_dict = generate_genomefile_dict(
        genomefile_B_f, int(args["--gapB"]), colour=True
    )

    # each alignment has coordinates to be recorded
    # store liftovers in a workdir for each pair
    alignment_coords = []
    alignment_coords = generate_alignment_dicts(
        liftover_f, genomefile_A_dict, genomefile_B_dict, alignment_coords
    )

    # current plot initiatior
    if i == 0:
        plot_chromosomes(
            genomefile_A_dict,
            1 - (2 * i),
            args["--labels"],
            int(args["--chromosome_width"]),
        )

    # for each alignment, if the two entries are from the two genomes, plot a curve between them
    for alignment in alignment_coords:
        # code is adapted from https://gist.github.com/andrewgiessel/5684769
        upper_x = alignment[0]
        lower_x = alignment[1]
        average_width = alignment[2]
        alignment_colour = alignment[3]
        Y_base = np.array(
            [-1, -0.568, -0.32, -0.16, -0.056, 0, 0.056, 0.16, 0.32, 0.568, 1]
        )
        Y = np.array([y - (2 * i) for y in Y_base])
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

        if args["--scaled"]:
            linewidth = linewidth_from_data_units(average_width, ax, reference="x")
        else:
            linewidth = float(args["--linewidth"])

        plt.plot(
            xnew,
            power_smooth,
            color=alignment_colour,
            alpha=float(args["--alpha"]),
            linewidth=linewidth,
        )

    # plot the chromosomes
    plot_chromosomes(
        genomefile_A_dict,
        1 - (2 * i),
        args["--labels"],
        int(args["--chromosome_width"]),
    )

    plot_chromosomes(
        genomefile_B_dict,
        -1 - (2 * i),
        args["--labels"],
        int(args["--chromosome_width"]),
    )


if __name__ == "__main__":
    args = docopt(__doc__)

    # read args and generate list of genome and busco files
    with open(args["--genomefiles"]) as fh:
        genomefiles = [line.rstrip() for line in fh]

    # set up plot
    fig = plt.figure(figsize=(28, 4 * (len(genomefiles) / 2)), frameon=False)
    ax = fig.add_subplot(111)
    ax.axis("off")
    plt.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis="y", which="both", right=False, left=False, labelleft=False)

    if args["--buscofiles"]:
        with open(args["--buscofiles"]) as fh:
            buscofiles = [line.rstrip() for line in fh]

        for i in range(len(genomefiles) - 1):
            plot_pair(
                i, genomefiles[i], genomefiles[i + 1], buscofiles[i], buscofiles[i + 1]
            )

    else:
        with open(args["--liftoverfiles"]) as fh:
            liftoverfiles = [line.rstrip() for line in fh]

        for i in range(len(genomefiles) - 1):
            plot_pair(i, genomefiles[i], genomefiles[i + 1], liftoverfiles[i])

    plt.savefig("busco3synteny.pdf", format="pdf", bbox_inches="tight")
    fig.set_frameon(True)
    plt.savefig(
        "busco3synteny.png",
        format="png",
        bbox_inches="tight",
        dpi=300,
        transparent=False,
        facecolor="w",
    )


"""
# example command
python busco3synteny.py -g genomefiles_2n.txt -b buscofiles2n.txt -t 0.25 -w 1
"""
