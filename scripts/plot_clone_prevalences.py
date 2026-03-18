from matplotlib import gridspec
from seaborn import color_palette

import colorcet
import json
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


def main(args):
    df = pd.read_table(args.in_df_file)

    colors = create_clone_color_dict(df.loc[~df["clone_id"].str.startswith("ancestral_"), "clone_id"].unique())

    with open(args.in_tree_file) as f:
        tree = nx.node_link_graph(json.load(f), edges="edges")

    setup_plot()

    fig = plt.figure(figsize=(10, 10), layout="constrained")

    gs = gridspec.GridSpec(1, 2, figure=fig)

    error_bars_ax = fig.add_subplot(gs[0, 0])

    plot_errorbars(colors, df, error_bars_ax)

    tree_ax = fig.add_subplot(gs[0, 1])

    plot_clone_tree(colors, tree, tree_ax)

    # fig.suptitle("Sample: {}".format(list(df["sample_id"].unique())[0]))

    fig.savefig(args.out_file, dpi=300, bbox_inches="tight")


def plot_errorbars(colors, df, error_bars_ax):
    df = df.loc[~df["clone_id"].str.startswith("ancestral_")]

    df = df.sort_values("clone_tree_order_idx", ascending=False)

    df["clone_color"] = df["clone_id"].map(colors)

    error_bars_ax = setup_axes(error_bars_ax)

    for i in df.index:
        x_val = df.loc[i, "mean_prevalence"]

        xerr_minus = x_val - df.loc[i, "lower_hdi"]

        xerr_minus = max(xerr_minus, 1e-6)

        xerr_plus = df.loc[i, "upper_hdi"] - x_val

        xerr_plus = max(xerr_plus, 1e-6)

        xerr = [[xerr_minus], [xerr_plus]]

        error_bars_ax.errorbar(
            y=df.loc[i, "clone_id"],
            x=df.loc[i, "mean_prevalence"],
            xerr=xerr,
            fmt='o',
            c=df.loc[i, "clone_color"],
        )

    error_bars_ax.set_xlabel("Clonal Prevalence")

    error_bars_ax.set_ylabel("Clone ID")


def plot_clone_tree(colors, tree, tree_ax):
    tree_ax.axis('off')

    pos = nx.nx_agraph.pygraphviz_layout(tree, prog="dot")

    node_border_colours = []

    node_colors = []

    node_sizes = []

    for n in tree.nodes:
        if n in colors:
            node_border_colours.append(colors[n])

            node_colors.append(colors[n])

        else:
            node_border_colours.append("k")

            node_colors.append("k")

        w = tree.nodes[n]["prevalence_stats"]["mean_prevalence"]

        node_sizes.append(512 * (w ** (3 / 2)))

    nx.draw_networkx(
        tree,
        pos=pos,
        ax=tree_ax,
        with_labels=False,
        node_size=node_sizes,
        node_color=node_colors,
        linewidths=0.75,
        edgecolors=node_border_colours,
    )


def create_clone_color_dict(node_names):
    num_colours = len(node_names)

    if num_colours > 10:
        colour_list = color_palette(colorcet.glasbey, num_colours).as_hex()

    else:
        colour_list = color_palette('colorblind', num_colours).as_hex()

    gt_colour_dict = dict(zip(node_names, colour_list))

    return gt_colour_dict


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

    parser.add_argument("-d", "--in-df-file", required=True)

    parser.add_argument("-t", "--in-tree-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
