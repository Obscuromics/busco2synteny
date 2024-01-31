#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: busco2synteny.py -a <STR> -b <STR> -x <STR> -y <STR> [-l <STR> -m <INT> -n <INT> -c <INT> -t <FLT> (-w <FLT> | -s) -z -h]

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
    -w, --line_width <FLT>                       Linewidth of alignments [default: 1]
    -s, --scaled                                Scale linewidths by aligned region length
    -z, --busco_colours                          Colour by BUSCO statuses (default: colour by ref chrom)
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
import functools
from string import ascii_lowercase

# mamba install docopt scipy numpy pandas matplotlib pyarrow

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
    i, liftover, genomefile_A_dict, genomefile_B_dict, alignment_coords
):
    colnums = [0] + [n + (3 * i) for n in range(1, 7)]
    liftover_df = pd.read_csv(
        liftover, sep="\t", usecols=colnums, header=None, index_col=None
    )
    # for each alignment, record coordinates in either genome
    for row in liftover_df.iterrows():
        alignment = []
        try:
            seqanc, seqA, startA, endA, seqB, startB, endB = row[1]
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
        except ValueError:
            next

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
                line_colour = cm.tab10((int(float(seqanc))-1))

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
    chromosomes = pd.read_csv(genomefile, sep='\t', usecols=[0], header=None)[0].to_list()
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
    labels = []
    for i in range(len(seqs)):
        labels.append(int(i + 1))
    return labels


def label_colours_by_ref(df):
    seqs = sorted(df["seq_a"].dropna().unique(), reverse=True)
    labels = get_labels(seqs)
    for seq, label in zip(seqs, labels):
        df.loc[df["seq_a"] == seq, "colour"] = label


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

    cols = ["colour"] + [
        label
        for label in df.columns
        if any(x in label for x in ["seq", "start", "stop"])
    ]

    # DROPNA LOCATION
    df.dropna(inplace=True)
    df[cols].to_csv("liftover.tsv", sep="\t", index=False, header=False)


if __name__ == "__main__":
    args = docopt(__doc__)

    # generate liftover file
    create_liftover_from_busco([args["--buscoA"], args["--buscoB"]],[args["--genomefileA"], args["--genomefileB"]])

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
    alignment_coords = generate_alignment_dicts(0, 
        "liftover.tsv", genomefile_A_dict, genomefile_B_dict, alignment_coords
    )

    # current plot initiatior
    plot_chromosomes(
        genomefile_A_dict,
        1,
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

        if args["--scaled"]:
            linewidth = linewidth_from_data_units(average_width, ax, reference="x")
        else:
            linewidth = float(args["--line_width"])

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
        genomefile_A_dict, 1, args["--labels"], int(args["--chromosome_width"])
    )
    plot_chromosomes(
        genomefile_B_dict, -1, args["--labels"], int(args["--chromosome_width"])
    )

    plt.savefig("busco2synteny.pdf", format="pdf", bbox_inches="tight")
    fig.set_frameon(True)
    plt.savefig(
        "busco2synteny.png",
        format="png",
        bbox_inches="tight",
        dpi=300,
        transparent=False,
        facecolor="w",
    )


"""
# example command
busco2synteny.py -a fol_ang.genomefile.tsv -b smi_aqu.genomefile.tsv -x busco_results/fol_ang/run_arthropoda_odb10/full_table.tsv -y busco_results/smi_aqu/run_arthropoda_odb10/full_table.tsv -t 0.25
"""
