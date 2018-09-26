#/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_args():
    parser = argparse.ArgumentParser(
        description="Wrapper of C++ bam_stats, add identity distplot",
        usage="usage: %(prog)s [options]")
    parser.add_argument("-b", "--bam",
        help="bam file [default %(default)s]", metavar="FILE")
    parser.add_argument("-p", "--prefix",
        help="output prefix [default %(default)s]", metavar="STR")
    parser.add_argument("-q", "--get_qual", default=False,
        help="get sequencing quality of the segment [default %(default)s]",
        metavar="FILE")
    if len(sys.argv) <= 1:
        parser.print_help()
        exit()
    return parser.parse_args()


def identity_plot(fragment_stats, prefix):
    stats = pd.read_csv(fragment_stats, sep="\t")
    pid = stats.FRAGMENT_IDENTITY
    num_bins = 50

    fig, ax = plt.subplots()

    # the histogram of the data
    n, bins, patches = ax.hist(pid, num_bins, density=1, edgecolor="black")

    ax.set_xlabel('Percent of identity')
    ax.set_ylabel('Probability density')
    # ax.set_title(r'Histogram of IQ: $\mu=100$, $\sigma=15$')

    # Tweak spacing to prevent clipping of ylabel
    fig.tight_layout()
    fig.set_size_inches(8,6)
    fig.savefig(prefix+".mapping.identity.png", dpi=300)

def main():
    args = get_args()
    cpp_bam_stats = os.path.join(os.path.dirname(os.path.abspath(__file__)),
            "bam_stats")
    subprocess.run([cpp_bam_stats, "-b", args.bam, "-p", args.prefix])
    fragment_stats_table = args.prefix+".mapping.fragment.stats.txt"
    identity_plot(fragment_stats_table, args.prefix)

if __name__ == "__main__":
    main()
