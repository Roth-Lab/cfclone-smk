import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import LinearLocator


colour_blind_friendly = {
    'blue': '#377eb8',
    'orange': '#ff7f00',
    'green': '#4daf4a',
    'pink': '#f781bf',
    'brown': '#a65628',
    'purple': '#984ea3',
    'grey': '#999999',
    'red': '#e41a1c',
    'yellow': '#dede00'
}


plot_settings = {
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "axes.titleweight": "bold",
    "axes.facecolor": "white",
    "grid.color": "darkgrey",
    "font.size": 10,
    "legend.title_fontsize": 16,
    "legend.fontsize": 12,
    "font.family": "sans-serif",
    "font.sans-serif": [
        'Helvetica',
        'Nimbus Sans',
        'Liberation Sans',
        'DejaVu Sans',
        'Arial',
    ],
}
    

def plot_samples(
    df: pd.DataFrame,
    figsize: tuple[int, int] | None = None,
    show: bool = False,
    output_path: str | None = None,
    pdf: PdfPages | None = None,
) -> None:
    
    with plt.rc_context(plot_settings):
    
        chains = df["chain"].unique().tolist()
        
        if len(chains) > len(colour_blind_friendly):
            
            raise ValueError("Number of chains exceeds number of available colours.")
        
        chain_colours = dict(zip(chains, colour_blind_friendly.values()))
        
        parameters = df["parameter"].unique().tolist()
        
        num_parameters = len(parameters)
        
        num_rows = num_parameters * 2
        
        if figsize is None:
            
            figsize = (15, num_rows * 2)

        fig = plt.figure(figsize=figsize, constrained_layout=True)
        
        gs = fig.add_gridspec(nrows=num_rows, ncols=2, height_ratios=[0.1, 0.9] * num_parameters)

        axes = np.empty(shape=(num_rows, 2), dtype=object)
        
        for param_idx in range(num_parameters):
            
            row_idx = 2 * param_idx
            
            axes[row_idx, 0] = fig.add_subplot(gs[row_idx, :])
            
            axes[row_idx + 1, 0] = fig.add_subplot(gs[row_idx + 1, 0])
            
            axes[row_idx + 1, 1] = fig.add_subplot(gs[row_idx + 1, 1])

        for param_idx in range(num_parameters):
            
            df_param = df[df["parameter"] == parameters[param_idx]]
            
            row_idx = 2 * param_idx

            axes[row_idx, 0].set_title("{}".format(parameters[param_idx]))
            
            axes[row_idx, 0].axis("off")

            for c in df_param["chain"].unique():
                
                df_c = df_param[df_param["chain"] == c]
                
                axes[row_idx + 1, 0].plot(
                    df_c["iteration"], 
                    df_c["sample"],
                    label=c,
                    alpha=0.25,
                    color=chain_colours[c]
                )
                
                axes[row_idx + 1, 1].hist(
                    df_c["sample"].values, 
                    bins=100, 
                    label=c, 
                    alpha=0.25, 
                    density=True, 
                    color=chain_colours[c]
                )

            axes[row_idx + 1, 0].spines["top"].set_visible(False)
            
            axes[row_idx + 1, 0].spines["right"].set_visible(False)
            
            axes[row_idx + 1, 0].set_ylabel("Sample value")
            
            axes[row_idx + 1, 0].set_xlabel("Iteration")

            axes[row_idx + 1, 1].spines["top"].set_visible(False)
            
            axes[row_idx + 1, 1].spines["right"].set_visible(False)
            
            axes[row_idx + 1, 1].set_ylabel("Density")
            
            axes[row_idx + 1, 1].set_xlabel("Sample value")

            axes[row_idx + 1, 0].legend(title="chain")
            
            axes[row_idx + 1, 1].legend(title="chain")

            axes[row_idx + 1, 0].yaxis.set_major_locator(LinearLocator(numticks=5))
            
            axes[row_idx + 1, 1].yaxis.set_major_locator(LinearLocator(numticks=5))

        fig.align_labels()

        if output_path is not None:
            
            plt.savefig(output_path, dpi=300)

        if pdf is not None:
            
            pdf.savefig(fig)

        if show:
            
            plt.show()

        plt.close()



def main(args):
    
    df_samples = (
        
        pd.read_csv(args.in_file, sep="\t")
        
        .melt(id_vars=["iteration", "chain"], var_name="parameter", value_name="sample")
        
        .filter(items=["parameter", "chain", "iteration", "sample"])
        
    )
    
    with PdfPages(args.out_file) as pdf:
        
        for param in df_samples["parameter"].unique():
            
            df_param = df_samples[df_samples["parameter"] == param]
            
            plot_samples(df=df_param, pdf=pdf)
            
            
    
        

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", type=str, required=True)

    parser.add_argument("-o", "--out-file", type=str, required=True)

    cli_args = parser.parse_args()

    main(cli_args)
