#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: busco3synteny.py -a <STR> (-x <STR> | -y <STR>)  [-l <STR> -m <INT> -c <INT> -t <FLT> (-w <FLT> | -s) -z -h]

  [Options]
    -a, --genomefiles <STR>                     File of genome file paths
    -x, --buscofiles <STR>                      File of busco file paths
    -y, --liftoverfile <STR>                    Liftover file
    -l, --labels <STR>                          Whether to plot labels, choose from True or False [default: True]
    -m, --gap <STR>                             Chromosome gap ratio per genome in units of 10Mb (eg '-m 1,1,1,1' for 4 genomes, default: equal)
    -c, --chromosome_width <INT>                Chromosome width [default: 6]
    -t, --alpha <FLT>                           Alpha of alignments [default: 0.1]
    -w, --linewidth <FLT>                       Linewidth of alignments [default: 1]
    -s, --scaled                                Scale linewidths by aligned region length
    -z, --busco_colours                         Colour by BUSCO statuses (default: colour by ref chrom) (not generalised in this script)
    -h, --help                                  Show this message

"""

# Example command:
# python scripts/busco3synteny.py -a genomefile_paths.txt -x arthropod.buscofile_paths.txt

import sys
import collections
import functools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from operator import itemgetter
from scipy.interpolate import make_interp_spline, BSpline
from string import ascii_lowercase
from docopt import docopt


# mamba install docopt scipy numpy pandas matplotlib pyarrow

# Generalising busco2synteny.py for N genomes
# Modified version of a script originally created by Alex Mackintosh
# https://github.com/A-J-F-Mackintosh/Mackintosh_et_al_2022_Binodaphne/blob/main/minimap2synteny.py


"""
Currently labels by one reference
    Ignore alignments between non-reference genomes
Clear nametags for intermediate genomes?
Add proper chromosome filtering by genomefile step
    Import the column and make a list of lists, check its there before merging
    Be good to add this to busco2synteny too
    Purely for colour labelling
    Maybe at the stage where it actually loads the liftover could do this?
"""


def generate_genomefile_dict(genomefile, offset, colour):
    genomefile_dict = {}
    cumulative_genome = 10_000_000 * offset
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
            cumulative_genome += 10_000_000 * offset
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
        if y_coord < -1:
            y = y_coord * 1.1
        else:
            y = y_coord * 1.2
        if labels == "True":
            plt.text(
                middle_of_chromosome,
                y,
                genomefile_dict[chromosome][3],
                ha="center",
                va="center",
                wrap=True,
                fontsize=20,
            )


def generate_alignment_dicts(
    i, liftover, genomefile_A_dict, genomefile_B_dict, alignment_coords
):
    colnums = [0] + [n + (3 * i) for n in range(1, 7)]
    liftover_df = pd.read_csv(
        liftover, sep="\t", usecols=colnums, header=None, index_col=None
    )
    # for each alignment, record coordinates in either genome
    for row in liftover_df.iterrows():
        alignment = []
        seqanc, seqA, startA, endA, seqB, startB, endB = row[1]
        if (seqA != seqA) or (seqB != seqB):
            continue
        midpointA = int(float(startA)) + (
            (int(float(endA)) - int(float(startA))) / 2
        )
        midpointB = int(float(startB)) + (
            (int(float(endB)) - int(float(startB))) / 2
        )
        average_width = (
            (int(float(endA)) - int(float(startA)))
            + (int(float(endB)) - int(float(startB)))
        ) / 2

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
            if int(float(seqanc)) == 0:
                line_colour = "lightgrey"
            else:
                line_colour = cm.tab20((int(float(seqanc)) - 1) / 9)

        # only interested in alignments on sequences in both genomefiles
        if len(alignment) == 2:
            alignment.append(average_width)
            alignment.append(line_colour)
            alignment.append(seqanc)
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


def load_busco_results(buscofile, genomefile):
    busco_colnames = ["busco_id", "status", "seq_code", "start", "stop"]
    chromosomes = pd.read_csv(genomefile, sep="\t", usecols=[0], header=None)[
        0
    ].to_list()
    df = pd.read_csv(
        buscofile,
        sep="\t",
        usecols=[0, 1, 2, 3, 4],
        header=None,
        names=busco_colnames,
        skiprows=[0, 1, 2],
    )

    df = df[df["status"] != "Missing"]
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
    df = df[df["seq"].isin(chromosomes)]
    return df[["busco_id", "status", "seq", "start", "stop"]]


def get_labels(seqs):
    labels = {}
    for i, seq in enumerate(seqs):
        labels[seq] = int(i + 1)
    labels[np.nan] = np.nan
    return labels


def label_colours_by_ref(df):
    seqs = []
    for genome in df.filter(like="seq").columns:
        seqs = seqs + sorted(df[genome].dropna().unique())
    labels = get_labels(seqs)
    df["colour"] = np.nan
    for i, genome in enumerate(df.filter(like="seq").columns):
        if (i + 1) == len(df.filter(like="seq").columns):
            break
        for seq in df[genome].dropna().unique():
            df.loc[
                (df[genome] == seq) & (df["colour"] != df["colour"]), "colour"
            ] = labels[seq]


def label_colours_by_busco(df):
    df["combs"] = df.filter(regex="status").astype(str).agg("_".join, axis=1)
    combs = np.sort(df["combs"].unique())
    labels = [n for n in range(len(combs))]
    global colour_dict
    colour_dict = {}
    for comb, label in zip(combs, labels):
        colour_dict[label] = comb
        df.loc[df["combs"] == comb, "colour"] = label
    df.drop(labels=["combs"], axis=1, inplace=True)


def create_liftover_from_busco(buscofile_list, genomefile_list):
    results_dfs = []
    for i, file in enumerate(buscofile_list):
        results_dfs.append(load_busco_results(file, genomefile_list[i]))

    it = iter(ascii_lowercase)

    for i, df in enumerate(results_dfs, start=1):
        df_label = next(it)
        df.rename(
            columns={
                col: "{}_{}".format(col, df_label)
                for col in ("seq", "start", "stop", "status")
            },
            inplace=True,
        )
    merge = functools.partial(pd.merge, how="outer", on="busco_id")
    df = functools.reduce(merge, results_dfs)

    if args["--busco_colours"]:
        label_colours_by_busco(df)
    else:
        label_colours_by_ref(df)

    # label_colours_by_ref(df)

    cols = ["colour"] + [
        label
        for label in df.columns
        if any(x in label for x in ["seq", "start", "stop"])
    ]

    df[cols].to_csv("liftover.tsv", sep="\t", index=False, header=False, na_rep="NA")


def plot_pair(
    i,
    genomefile_A_f,
    genomefile_B_f,
    liftover_f="liftover.tsv",
):
    # generate dicts for each genome with cumulative coordinates
    if args["--gap"]:
        gapA = int(gap_ratios[i])
        gapB = int(gap_ratios[i + 1])
    else:
        gapA = gapB = 1

    genomefile_A_dict = generate_genomefile_dict(genomefile_A_f, gapA, colour=False)
    genomefile_B_dict = generate_genomefile_dict(genomefile_B_f, gapB, colour=True)

    # each alignment has coordinates to be recorded
    # store liftovers in a workdir for each pair
    alignment_coords = []
    alignment_coords = generate_alignment_dicts(
        i, liftover_f, genomefile_A_dict, genomefile_B_dict, alignment_coords
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
    labelled = []
    for alignment in alignment_coords:
        # code is adapted from https://gist.github.com/andrewgiessel/5684769
        upper_x = alignment[0]
        lower_x = alignment[1]
        average_width = alignment[2]
        alignment_colour = alignment[3]
        colour_code = alignment[4]
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

        if args["--busco_colours"]:
            if alignment_colour not in labelled:
                labelled.append(alignment_colour)

                plt.plot(
                    xnew,
                    power_smooth,
                    color=alignment_colour,
                    label=colour_dict[float(colour_code)],
                    alpha=float(args["--alpha"]),
                    linewidth=linewidth,
                )

            else:
                plt.plot(
                    xnew,
                    power_smooth,
                    color=alignment_colour,
                    alpha=float(args["--alpha"]),
                    linewidth=linewidth,
                )
            plt.legend()

        else:
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

    if args["--gap"]:
        gap_ratios = str(args["--gap"]).split(",")

    if args["--buscofiles"]:
        with open(args["--buscofiles"]) as fh:
            buscofiles = [line.rstrip() for line in fh]

        create_liftover_from_busco(buscofiles, genomefiles)

        for i in range(len(genomefiles) - 1):
            plot_pair(i, genomefiles[i], genomefiles[i + 1])

    else:
        for i in range(len(genomefiles) - 1):
            plot_pair(i, genomefiles[i], genomefiles[i + 1], args["--liftoverfile"])

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
