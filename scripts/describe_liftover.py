#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: describe_liftover.py -l <STR> [-h]

  [Options]
    -l, --liftover <STR>                     Liftover file

"""

# describe_liftover.py -f liftover.tsv

from docopt import docopt
import pandas as pd
import numpy as np
from string import ascii_lowercase
import itertools
from scipy.stats import spearmanr

pd.options.mode.chained_assignment = None  # default='warn'


def get_genome_df(df, genome):
    gen_loc = df.columns.get_loc(genome)
    genome_df = df.iloc[:, gen_loc : gen_loc + 3]
    return genome_df


def compare_chroms(df, a_chr, b_chr):
    tmp_df = df[(df.iloc[:, 0] == a_chr) & (df.iloc[:, 3] == b_chr)]
    chr_df = tmp_df.drop_duplicates()

    a_aln_len = int(np.sum(chr_df.iloc[:, 2] - chr_df.iloc[:, 1]))
    b_aln_len = int(np.sum(chr_df.iloc[:, 5] - chr_df.iloc[:, 4]))
    a_prop = a_aln_len / length_dict[a_chr]
    b_prop = b_aln_len / length_dict[b_chr]

    chr_df["a_mid"] = chr_df.iloc[:, 1] + ((chr_df.iloc[:, 2] - chr_df.iloc[:, 1]) / 2)
    chr_df["b_mid"] = chr_df.iloc[:, 4] + ((chr_df.iloc[:, 5] - chr_df.iloc[:, 4]) / 2)
    chr_df = chr_df.sort_values(by="a_mid")

    b_by_a = chr_df["b_mid"]
    b_by_b = np.sort(b_by_a)
    corr = np.abs(spearmanr(b_by_a, b_by_b).statistic)
    return a_prop, b_prop, corr


def compare_genomes(genome_a, genome_b):
    a_df = get_genome_df(df, genome_a)
    b_df = get_genome_df(df, genome_b)

    a_seqs = a_df.iloc[:, 0].dropna().unique()
    b_seqs = b_df.iloc[:, 0].dropna().unique()

    global length_dict
    length_dict = {}
    for seq in a_seqs:
        length_dict[seq] = int(
            np.sum(
                a_df[a_df.iloc[:, 0] == seq].iloc[:, 2]
                - a_df[a_df.iloc[:, 0] == seq].iloc[:, 1]
            )
        )
    for seq in b_seqs:
        length_dict[seq] = int(
            np.sum(
                b_df[b_df.iloc[:, 0] == seq].iloc[:, 2]
                - b_df[b_df.iloc[:, 0] == seq].iloc[:, 1]
            )
        )

    comp_df = a_df.join(b_df)

    combs = list(itertools.product(a_seqs, b_seqs))
    test_comb = combs[0]

    results = []
    for a_chr, b_chr in combs:
        a_prop, b_prop, corr = compare_chroms(comp_df, a_chr, b_chr)
        results.append([a_chr, b_chr, a_prop, b_prop, corr])

    columns = ["chr_a", "chr_b", "prop_a", "prop_b", "corr"]
    results_df = pd.DataFrame(results, columns=columns)
    results_df["genome_a"] = [genome_a] * results_df.shape[0]
    results_df["genome_b"] = [genome_b] * results_df.shape[0]

    return results_df[["genome_a", "genome_b"] + columns]


if __name__ == "__main__":
    args = docopt(__doc__)
    liftover_f = args["--liftover"]

    df = pd.read_csv(liftover_f, sep="\t", header=None)
    df.drop(0, axis=1, inplace=True)
    n = int(len(df.columns) / 3)

    it = iter(ascii_lowercase)
    colnames = []
    for i in range(n):
        label = next(it)
        colnames += [prefix % label for prefix in ["seq_%s", "start_%s", "stop_%s"]]
    df.columns = colnames

    genomes = df.filter(like="seq").columns.to_list()
    gen_combs = list(itertools.combinations(genomes, 2))

    out_dfs = []
    for gen_a, gen_b in gen_combs:
        out_dfs.append(compare_genomes(gen_a, gen_b))

    out_df = pd.concat(out_dfs)
    out_df.reset_index().drop("index", axis=1).to_csv(
        "liftover_stats.tsv", sep="\t", index=False
    )
