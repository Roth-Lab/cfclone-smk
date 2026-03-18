from matplotlib.collections import LineCollection
from matplotlib import gridspec
from matplotlib.patches import Patch
from seaborn.matrix import ClusterGrid
from seaborn.utils import despine, to_utf8, axis_ticklabels_overlap
from skbio.tree import TreeNode

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings


def main(args):
    df = pd.read_csv(args.in_file, index_col="clone_id", sep="\t")

    df.index = df.index.astype(str)

    for i in range(df.shape[0]):
        df.iloc[i, i] = np.nan

    tree = TreeNode.read(args.tree_file, convert_underscores=False)

    cmap_greater_than = mcolors.LinearSegmentedColormap.from_list(
        "red_white_blue", [(0, "blue"), (0.5, "white"), (1, "red")]
    )

    fig = clustermap_phylo(
        df,
        cmap=cmap_greater_than,
        cbar_kws={'label': r"$P(\rho_{i} > \rho_{j})$"},
        col_tree=tree,
        row_tree=tree,
        use_branch_lengths=False,
        label=True,
        vmin=0,
        vmax=1,
    )

    fig.ax_heatmap.set(xlabel='Clone ID: $\\mathbf{i}$', ylabel='Clone ID: $\\mathbf{j}$')

    if args.sample_id is not None:
        fig.ax_col_dendrogram.set(title="Sample: {}".format(args.sample_id))

    plt.tight_layout()

    plt.savefig(args.out_file, dpi=300, bbox_inches='tight')

    plt.close()


def _draw_figure(fig):
    """Force draw of a matplotlib figure, accounting for back-compat."""
    # See https://github.com/matplotlib/matplotlib/issues/19197 for context
    fig.canvas.draw()
    if fig.stale:
        try:
            fig.draw(fig.canvas.get_renderer())
        except AttributeError:
            pass


def _index_to_label(index):
    """Convert a pandas index or multiindex to an axis label."""
    if isinstance(index, pd.MultiIndex):
        return "-".join(map(to_utf8, index.names))
    else:
        return index.name


def _index_to_ticklabels(index):
    """Convert a pandas index or multiindex into ticklabels."""
    if isinstance(index, pd.MultiIndex):
        return ["-".join(map(to_utf8, i)) for i in index.values]
    else:
        return index.values


class DendrogramPhyloPlotter:
    """Object for drawing tree of similarities between data rows/columns"""

    def __init__(self, data, tree, metric, method, axis, label, rotate, use_branch_lengths):
        """Plot a dendrogram of the relationships between the columns of data

        Parameters
        ----------
        data : pandas.DataFrame
            Rectangular data
        """
        self.axis = axis
        if self.axis == 1:
            data = data.T

        if isinstance(data, pd.DataFrame):
            array = data.values
        else:
            array = np.asarray(data)
            data = pd.DataFrame(array)

        self.array = array
        self.data = data

        self.shape = self.data.shape
        self.metric = metric
        self.method = method
        self.axis = axis
        self.label = label
        self.rotate = rotate
        self.tree = tree
        self.use_branch_lengths = use_branch_lengths

        self._preprocess_tree_to_extant_tips()

        self.dendrogram = self.calculate_dendrogram()

        ticks = np.arange(self.data.shape[0])

        if self.label:
            ticklabels = _index_to_ticklabels(self.data.index)
            ticklabels = [ticklabels[i] for i in self.reordered_ind]
            if self.rotate:
                self.xticks = []
                self.yticks = ticks
                self.xticklabels = []

                self.yticklabels = ticklabels
                self.ylabel = _index_to_label(self.data.index)
                self.xlabel = ''
            else:
                self.xticks = ticks
                self.yticks = []
                self.xticklabels = ticklabels
                self.yticklabels = []
                self.ylabel = ''
                self.xlabel = _index_to_label(self.data.index)
        else:
            self.xticks, self.yticks = [], []
            self.yticklabels, self.xticklabels = [], []
            self.xlabel, self.ylabel = '', ''

        self.dependent_coord = self.dendrogram['dcoord']
        self.independent_coord = self.dendrogram['icoord']

    def calculate_dendrogram(self):
        """Calculates a dendrogram based on the linkage matrix

        Made a separate function, not a property because don't want to
        recalculate the dendrogram every time it is accessed.

        Returns
        -------
        dendrogram : dict
            Dendrogram dictionary as returned by scipy.cluster.hierarchy
            .dendrogram. The important key-value pairing is
            "reordered_ind" which indicates the re-ordering of the matrix
        """

        dendrogram = dict()

        use_length = self.use_branch_lengths

        icoords = []
        dcoords = []

        node_coord_dict = {}

        x_curr = 0.5

        tree = self.tree

        x_spacing = 1

        for node in tree.postorder():
            curr_height, _ = node.height(include_self=True, missing_as_zero=True, use_length=use_length)

            node_y = curr_height

            if node.is_tip():
                node_x = x_curr * x_spacing
                x_curr += 1

                icoords.append([node_x, node_x])
                dcoords.append([node_y, 0])

            else:
                child_xcoords = np.array([node_coord_dict[child][0] for child in node.children])
                node_x = child_xcoords.mean()

                child_ycoords = [node_coord_dict[child][1] for child in node.children]

                icoords.append([child_xcoords.min(), child_xcoords.max()])
                dcoords.append([node_y, node_y])

                for child_x, child_y in zip(child_xcoords, child_ycoords):
                    icoords.append([child_x, child_x])
                    dcoords.append([node_y, child_y])

            node_coord_dict[node] = (node_x, node_y)

        dendrogram["icoord"] = icoords
        dendrogram["dcoord"] = dcoords

        dendrogram_leaves = list(self.data.index.get_loc(tip.name) for tip in tree.tips())
        missing_idx = list(self.data.index.get_loc(missing) for missing in self.missing_indices)
        dendrogram_leaves.extend(missing_idx)
        dendrogram['leaves'] = np.array(dendrogram_leaves)

        return dendrogram

    def _preprocess_tree_to_extant_tips(self):
        index_set = set(self.data.index)
        tip_set = set(tip.name for tip in self.tree.tips())
        tip_intersect = tip_set.intersection(index_set)
        self.missing_indices = index_set - tip_set
        self.tree = self.tree.shear(tip_intersect, prune=True, strict=False, inplace=False)

    @property
    def reordered_ind(self):
        """Indices of the matrix, reordered by the dendrogram"""
        return self.dendrogram['leaves']

    def plot(self, ax, tree_kws):
        """Plots a dendrogram of the similarities between data on the axes

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object upon which the dendrogram is plotted

        """
        tree_kws = {} if tree_kws is None else tree_kws.copy()
        tree_kws.setdefault("linewidths", 0.5)
        tree_kws.setdefault("colors", tree_kws.pop("color", (0.2, 0.2, 0.2)))

        if self.rotate and self.axis == 0:
            coords = zip(self.dependent_coord, self.independent_coord)
        else:
            coords = zip(self.independent_coord, self.dependent_coord)
        lines = LineCollection([list(zip(x, y)) for x, y in coords], **tree_kws)

        ax.add_collection(lines)
        number_of_leaves = len(self.reordered_ind) - len(self.missing_indices)
        max_dependent_coord = max(map(max, self.dependent_coord))

        if self.rotate:
            ax.yaxis.set_ticks_position('right')
            ax.set_xlim(0, max_dependent_coord)
            ax.set_ylim(0, number_of_leaves)
            ax.invert_xaxis()
            ax.invert_yaxis()

        else:
            ax.set_xlim(0, number_of_leaves)
            ax.set_ylim(0, max_dependent_coord)

        despine(ax=ax, bottom=True, left=True)

        ax.set(
            xticks=self.xticks,
            yticks=self.yticks,
            xlabel=self.xlabel,
            ylabel=self.ylabel,
        )
        xtl = ax.set_xticklabels(self.xticklabels)
        ytl = ax.set_yticklabels(self.yticklabels, rotation='vertical')

        # Force a draw of the plot to avoid matplotlib window error
        _draw_figure(ax.figure)

        if len(ytl) > 0 and axis_ticklabels_overlap(ytl):
            plt.setp(ytl, rotation="horizontal")
        if len(xtl) > 0 and axis_ticklabels_overlap(xtl):
            plt.setp(xtl, rotation="vertical")
        return self


class ClusterGridPhylogeny(ClusterGrid):

    def __init__(
        self,
        data,
        pivot_kws=None,
        z_score=None,
        standard_scale=None,
        figsize=None,
        row_colors=None,
        col_colors=None,
        mask=None,
        dendrogram_ratio=None,
        colors_ratio=None,
        cbar_pos=None,
        legend_ratio=0.18,
        legend_luts=None,
    ):
        super().__init__(
            data,
            pivot_kws,
            z_score,
            standard_scale,
            figsize,
            row_colors,
            col_colors,
            mask,
            dendrogram_ratio,
            colors_ratio,
            cbar_pos,
        )

        self.is_legend_plotting_active = legend_luts is not None

        if self.is_legend_plotting_active:

            if legend_ratio == 0.0:
                legend_ratio = 0.1

            self.legend_luts = legend_luts

            plt.close(self._figure)

            self._figure = plt.figure(figsize=figsize)

            try:
                row_dendrogram_ratio, col_dendrogram_ratio = dendrogram_ratio
            except TypeError:
                row_dendrogram_ratio = col_dendrogram_ratio = dendrogram_ratio

            try:
                row_colors_ratio, col_colors_ratio = colors_ratio
            except TypeError:
                row_colors_ratio = col_colors_ratio = colors_ratio

            width_ratios = self.dim_ratios_legend(self.row_colors, row_dendrogram_ratio, row_colors_ratio, legend_ratio)

            height_ratios = self.dim_ratios(self.col_colors, col_dendrogram_ratio, col_colors_ratio)

            nrows = 2 if self.col_colors is None else 3
            ncols = 2 if self.row_colors is None else 3

            ncols += 1

            self.gs = gridspec.GridSpec(
                nrows,
                ncols,
                width_ratios=width_ratios,
                height_ratios=height_ratios,
                wspace=0,
            )

            self.ax_row_dendrogram = self._figure.add_subplot(self.gs[-1, 0])
            self.ax_col_dendrogram = self._figure.add_subplot(self.gs[0, -2])
            self.ax_row_dendrogram.set_axis_off()
            self.ax_col_dendrogram.set_axis_off()

            self.ax_row_colors = None
            self.ax_col_colors = None

            if self.row_colors is not None:
                self.ax_row_colors = self._figure.add_subplot(self.gs[-1, 1])
            if self.col_colors is not None:
                self.ax_col_colors = self._figure.add_subplot(self.gs[1, -2])

            self.ax_heatmap = self._figure.add_subplot(self.gs[-1, -2])

            if cbar_pos is None:
                self.ax_cbar = self.cax = None
            else:
                # Initialize the colorbar axes in the gridspec so that tight_layout
                # works. We will move it where it belongs later. This is a hack.
                self.ax_cbar = self._figure.add_subplot(self.gs[0, 0])
                self.cax = self.ax_cbar  # Backwards compatibility

            self.ax_legend = self._figure.add_subplot(self.gs[:, -1])
            self.ax_legend.set_axis_off()

    def dim_ratios_legend(self, colors, dendrogram_ratio, colors_ratio, legend_ratio):
        """Get the proportions of the figure taken up by each axes."""
        ratios = [dendrogram_ratio]

        if colors is not None:
            # Colors are encoded as rgb, so there is an extra dimension
            if np.ndim(colors) > 2:
                n_colors = len(colors)
            else:
                n_colors = 1

            ratios += [n_colors * colors_ratio]

        ratios_to_sum = ratios + [legend_ratio]
        ratios.append(1 - sum(ratios_to_sum))
        ratios.append(legend_ratio)

        return ratios

    def plot_dendrograms_phylo(
        self,
        row_cluster,
        col_cluster,
        metric,
        method,
        row_tree,
        col_tree,
        tree_kws,
        use_branch_lengths,
    ):
        if row_cluster:
            self.dendrogram_row = dendrogram_phylo(
                self.data2d,
                metric=metric,
                method=method,
                label=False,
                axis=0,
                ax=self.ax_row_dendrogram,
                rotate=True,
                tree=row_tree,
                tree_kws=tree_kws,
                use_branch_lengths=use_branch_lengths,
            )
        else:
            self.ax_row_dendrogram.set_xticks([])
            self.ax_row_dendrogram.set_yticks([])
        if col_cluster:
            self.dendrogram_col = dendrogram_phylo(
                self.data2d,
                metric=metric,
                method=method,
                label=False,
                axis=1,
                ax=self.ax_col_dendrogram,
                tree=col_tree,
                tree_kws=tree_kws,
                use_branch_lengths=use_branch_lengths,
            )
        else:
            self.ax_col_dendrogram.set_xticks([])
            self.ax_col_dendrogram.set_yticks([])
        despine(ax=self.ax_row_dendrogram, bottom=True, left=True)
        despine(ax=self.ax_col_dendrogram, bottom=True, left=True)

    def plot_phylo(
        self,
        metric,
        method,
        colorbar_kws,
        row_cluster,
        col_cluster,
        tree_kws,
        row_tree,
        col_tree,
        use_branch_lengths,
        **kws,
    ):

        # heatmap square=True sets the aspect ratio on the axes, but that is
        # not compatible with the multi-axes layout of clustergrid
        if kws.get("square", False):
            msg = "``square=True`` ignored in clustermap"
            warnings.warn(msg)
            kws.pop("square")

        colorbar_kws = {} if colorbar_kws is None else colorbar_kws

        self.plot_dendrograms_phylo(
            row_cluster,
            col_cluster,
            metric,
            method,
            row_tree=row_tree,
            col_tree=col_tree,
            tree_kws=tree_kws,
            use_branch_lengths=use_branch_lengths,
        )
        try:
            xind = self.dendrogram_col.reordered_ind
        except AttributeError:
            xind = np.arange(self.data2d.shape[1])
        try:
            yind = self.dendrogram_row.reordered_ind
        except AttributeError:
            yind = np.arange(self.data2d.shape[0])

        self.plot_colors(xind, yind, **kws)
        self.plot_matrix(colorbar_kws, xind, yind, **kws)
        if self.is_legend_plotting_active:
            self.plot_legends()
        return self

    def plot_legends(self):
        legend_ax = self.ax_legend
        for lut_dict in self.legend_luts:
            lut = lut_dict.pop('lut')
            lut_names = lut.keys()
            handles = [Patch(facecolor=lut[name]) for name in lut_names]
            curr_legend = legend_ax.legend(handles, lut_names, **lut_dict)
            legend_ax.add_artist(curr_legend)


def clustermap_phylo(
    data,
    *,
    pivot_kws=None,
    method='average',
    metric='euclidean',
    z_score=None,
    standard_scale=None,
    figsize=(10, 10),
    cbar_kws=None,
    row_cluster=True,
    col_cluster=True,
    row_colors=None,
    col_colors=None,
    mask=None,
    dendrogram_ratio=0.2,
    colors_ratio=0.03,
    cbar_pos=(0.02, 0.8, 0.05, 0.18),
    tree_kws=None,
    row_tree=None,
    col_tree=None,
    legend_ratio=0.18,
    legend_luts=None,
    use_branch_lengths=True,
    **kwargs,
):
    """
    Plot a matrix dataset as a hierarchically-clustered heatmap.

    This function requires scipy to be available.

    Parameters
    ----------
    data : 2D array-like
        Rectangular data for clustering. Cannot contain NAs.
    pivot_kws : dict, optional
        If `data` is a tidy dataframe, can provide keyword arguments for
        pivot to create a rectangular dataframe.
    method : str, optional
        Linkage method to use for calculating clusters. See
        :func:`scipy.cluster.hierarchy.linkage` documentation for more
        information.
    metric : str, optional
        Distance metric to use for the data. See
        :func:`scipy.spatial.distance.pdist` documentation for more options.
        To use different metrics (or methods) for rows and columns, you may
        construct each linkage matrix yourself and provide them as
        `{row,col}_linkage`.
    z_score : int or None, optional
        Either 0 (rows) or 1 (columns). Whether or not to calculate z-scores
        for the rows or the columns. Z scores are: z = (x - mean)/std, so
        values in each row (column) will get the mean of the row (column)
        subtracted, then divided by the standard deviation of the row (column).
        This ensures that each row (column) has mean of 0 and variance of 1.
    standard_scale : int or None, optional
        Either 0 (rows) or 1 (columns). Whether or not to standardize that
        dimension, meaning for each row or column, subtract the minimum and
        divide each by its maximum.
    figsize : tuple of (width, height), optional
        Overall size of the figure.
    cbar_kws : dict, optional
        Keyword arguments to pass to `cbar_kws` in :func:`heatmap`, e.g. to
        add a label to the colorbar.
    {row,col}_cluster : bool, optional
        If ``True``, cluster the {rows, columns}.
    {row,col}_linkage : :class:`numpy.ndarray`, optional
        Precomputed linkage matrix for the rows or columns. See
        :func:`scipy.cluster.hierarchy.linkage` for specific formats.
    {row,col}_colors : list-like or pandas DataFrame/Series, optional
        List of colors to label for either the rows or columns. Useful to evaluate
        whether samples within a group are clustered together. Can use nested lists or
        DataFrame for multiple color levels of labeling. If given as a
        :class:`pandas.DataFrame` or :class:`pandas.Series`, labels for the colors are
        extracted from the DataFrames column names or from the name of the Series.
        DataFrame/Series colors are also matched to the data by their index, ensuring
        colors are drawn in the correct order.
    mask : bool array or DataFrame, optional
        If passed, data will not be shown in cells where `mask` is True.
        Cells with missing values are automatically masked. Only used for
        visualizing, not for calculating.
    {dendrogram,colors}_ratio : float, or pair of floats, optional
        Proportion of the figure size devoted to the two marginal elements. If
        a pair is given, they correspond to (row, col) ratios.
    cbar_pos : tuple of (left, bottom, width, height), optional
        Position of the colorbar axes in the figure. Setting to ``None`` will
        disable the colorbar.
    tree_kws : dict, optional
        Parameters for the :class:`matplotlib.collections.LineCollection`
        that is used to plot the lines of the dendrogram tree.
    kwargs : other keyword arguments
        All other keyword arguments are passed to :func:`heatmap`.

    Returns
    -------
    :class:`ClusterGridPhylogeny`
        A :class:`ClusterGridPhylogeny` instance.

    See Also
    --------
    heatmap : Plot rectangular data as a color-encoded matrix.

    Notes
    -----
    The returned object has a ``savefig`` method that should be used if you
    want to save the figure object without clipping the dendrograms.

    To access the reordered row indices, use:
    ``clustergrid.dendrogram_row.reordered_ind``

    Column indices, use:
    ``clustergrid.dendrogram_col.reordered_ind``

    Examples
    --------

    .. include:: ../docstrings/clustermap.rst

    """
    plotter = ClusterGridPhylogeny(
        data,
        pivot_kws=pivot_kws,
        figsize=figsize,
        row_colors=row_colors,
        col_colors=col_colors,
        z_score=z_score,
        standard_scale=standard_scale,
        mask=mask,
        dendrogram_ratio=dendrogram_ratio,
        colors_ratio=colors_ratio,
        cbar_pos=cbar_pos,
        legend_ratio=legend_ratio,
        legend_luts=legend_luts,
    )

    return plotter.plot_phylo(
        metric=metric,
        method=method,
        colorbar_kws=cbar_kws,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        tree_kws=tree_kws,
        row_tree=row_tree,
        col_tree=col_tree,
        use_branch_lengths=use_branch_lengths,
        **kwargs,
    )


def dendrogram_phylo(
    data,
    tree,
    *,
    axis=1,
    label=True,
    metric='euclidean',
    method='average',
    rotate=False,
    tree_kws=None,
    ax=None,
    use_branch_lengths=True,
):
    """Draw a tree diagram of relationships within a matrix

    Parameters
    ----------
    data : pandas.DataFrame
        Rectangular data
    linkage : numpy.array, optional
        Linkage matrix
    axis : int, optional
        Which axis to use to calculate linkage. 0 is rows, 1 is columns.
    label : bool, optional
        If True, label the dendrogram at leaves with column or row names
    metric : str, optional
        Distance metric. Anything valid for scipy.spatial.distance.pdist
    method : str, optional
        Linkage method to use. Anything valid for
        scipy.cluster.hierarchy.linkage
    rotate : bool, optional
        When plotting the matrix, whether to rotate it 90 degrees
        counter-clockwise, so the leaves face right
    tree_kws : dict, optional
        Keyword arguments for the ``matplotlib.collections.LineCollection``
        that is used for plotting the lines of the dendrogram tree.
    ax : matplotlib axis, optional
        Axis to plot on, otherwise uses current axis

    Returns
    -------
    dendrogramplotter : _DendrogramPlotter
        A Dendrogram plotter object.

    Notes
    -----
    Access the reordered dendrogram indices with
    dendrogramplotter.reordered_ind

    """
    plotter = DendrogramPhyloPlotter(
        data,
        tree,
        axis=axis,
        metric=metric,
        method=method,
        label=label,
        rotate=rotate,
        use_branch_lengths=use_branch_lengths,
    )
    if ax is None:
        ax = plt.gca()

    return plotter.plot(ax=ax, tree_kws=tree_kws)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-t", "--tree-file", required=True)

    parser.add_argument("-s", "--sample-id", default=None)

    cli_args = parser.parse_args()

    main(cli_args)
