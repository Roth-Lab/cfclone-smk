from matplotlib.patches import Patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb


def main(args):
    setup_plot()

    fig = plt.figure(figsize=(16, 12))

    grid = fig.add_gridspec(6, 1, hspace=0.1)

    df = pd.read_table(args.in_file)

    chroms = sort_chroms(df["chrom"].unique())

    chroms_size = df["chrom"].value_counts()

    width_ratios = [chroms_size[x] for x in chroms]

    plot_vals = [
        ("data_rdr", "mu", "RDR", df, True),
        (None, "mu_residual", "RDR residual", df, False),
        (None, "rdr_outlier_prob", "RDR outlier probability", df, False),
        ("data_baf", "p", "BAF", df, False),
        (None, "p_residual", "BAF residual", df, False),
        (None, "baf_outlier_prob", "BAF outlier probability", df, False),
    ]
    for i, v in enumerate(plot_vals):
        sub_grid = grid[i].subgridspec(1, len(chroms), width_ratios=width_ratios, wspace=0.05)

        plot_fit(
            v[3],
            chroms,
            fig,
            sub_grid,
            title=v[2],
            y_col=v[0],
            y_col_fit=v[1],
            plot_legend=v[4],
        )

    # fig.suptitle("Sample: {}".format(sample_id))

    fig.supxlabel("Chromosome")

    grid.tight_layout(fig)

    fig.savefig(args.out_file, dpi=300, bbox_inches="tight")


def build_merged_param_summary_df(ctdna_df, param_summary_file):
    param_summary_df = pd.read_table(param_summary_file)
    param_summary_df = pd.merge(param_summary_df, ctdna_df, on="bin")
    return param_summary_df


def plot_fit(
    df,
    chroms,
    fig,
    grid,
    title=None,
    y_col="data_rdr",
    y_col_fit="mu",
    plot_legend=False,
):
    mean_col = y_col_fit + "_mean"

    lower_col = y_col_fit + "_lower_hdi"

    upper_col = y_col_fit + "_upper_hdi"

    y_max = df[upper_col].max()

    y_min = df[lower_col].min()

    if y_col is not None:
        y_max = max(df[y_col].max(), y_max)

        y_min = min(df[y_col].min(), y_min)

    weeyums_blue = "#1868DB"

    papaya_orng = "#F47600"

    grouped = df.groupby("chrom")

    last_chrom = chroms[-1]

    for i, chrom in enumerate(chroms):
        chrom_df = grouped.get_group(chrom)

        chrom_df = chrom_df.sort_values(by=["start"])

        num_bins = chrom_df.shape[0]

        chrom_df["idx"] = np.arange(num_bins)

        ax = fig.add_subplot(grid[0, i])

        setup_axes(ax)

        if y_col is not None:
            ax.scatter(
                chrom_df["idx"],
                chrom_df[y_col],
                c=papaya_orng,
                s=1,
                alpha=0.5,
            )

        ax.scatter(
            chrom_df["idx"],
            chrom_df[mean_col],
            c=weeyums_blue,
            s=1,
            alpha=0.5,
        )

        ax.fill_between(
            chrom_df["idx"],
            y1=chrom_df[lower_col],
            y2=chrom_df[upper_col],
            color=weeyums_blue,
            alpha=0.2,
        )

        sb.despine(ax=ax, offset=10)

        ax.spines["top"].set_visible(False)

        ax.spines["right"].set_visible(False)

        ax.xaxis.grid(False)

        if i != 0:
            ax.spines["left"].set_visible(False)

            ax.tick_params(axis='y', labelleft=False, left=False)

        else:
            ax.tick_params(axis="x", which="major", labelsize=12)

        ax.set_xticks([num_bins / 2])

        ax.set_xticklabels([chrom.replace("chr", "")], fontsize=12)

        ax.set_ylim(y_min, y_max)

        if plot_legend and chrom == last_chrom:
            handles = [Patch(facecolor=papaya_orng), Patch(facecolor=weeyums_blue)]

            labels = ("Data", "Model")

            ax.legend(handles, labels, loc="upper right", bbox_to_anchor=(1, 1.15))

    if title is not None:
        ax = fig.add_subplot(grid[:])

        ax.axis("off")

        ax.set_title(title)


# sort_chroms adapted from: https://github.com/Roth-Lab/hapclone-smk/blob/main/scripts/plot_clone_pseudobulk.py
def sort_chroms(chroms):
    numeric = []
    string = []

    if chroms[0].startswith("chr"):
        chr_prefix = True
    else:
        chr_prefix = False

    for c in chroms:
        if chr_prefix:
            c = c.replace("chr", "")
        try:
            numeric.append(int(c))
        except ValueError:
            string.append(c)

    chroms = [str(x) for x in sorted(numeric)] + list(sorted(string))

    if chr_prefix:
        chroms = ["chr{}".format(x) for x in chroms]

    return chroms


def setup_plot(font_size=None):
    """
    Set up matplotlib in a consistent way between plots.
    """
    plt.rcParams['font.family'] = ' sans-serif'
    plt.rcParams['font.sans-serif'] = [
        'Helvetica',
        'Nimbus Sans',
        'Liberation Sans',
        'DejaVu Sans',
        'Arial',
    ]
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['grid.color'] = 'darkgrey'
    if font_size is not None:
        plt.rcParams["font.size"] = font_size


def setup_axes(ax):
    """
    Set up axes in a consistent way between plots.  Spines are at left and bottom, offset by 10pts.
    Other spines not visible.  Ticks only for left and bottom.
    """
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()

    ax.xaxis.grid(True, which="major", linestyle=':')
    ax.yaxis.grid(True, which="major", linestyle=':')

    return ax


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
