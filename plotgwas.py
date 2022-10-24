import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import sys
import argparse

"""
Run examples:
(1) Manhattan plot with one sumstats:
python plotgwas.py manhattan --config config.plotgwas.3.cfg --out manhattan.3.svg
(2) Miami plot with two sumstats in the top panel and one sumstats in the bottom panel saved as png and as svg files:
python plotgwas.py miami --config-top config.plotgwas.1.cfg config.plotgwas.2.cfg --config-bottom config.plotgwas.3.cfg --out miami.1.2.3.png miami.1.2.3.svg
"""


def parse_args(args):
    parser = argparse.ArgumentParser(description="Tools to plot GWAS summary statistic data.")
    subparsers = parser.add_subparsers()

    parser_manhattan = subparsers.add_parser("manhattan", help="Make Manhattan plot.")
    parser_manhattan.add_argument("--config", type=str, nargs='+', required=True, help="List of config files to plot.")
    parser_manhattan.add_argument("--out", type=str, nargs='+', required=True, help="List of output file names.")
    parser_manhattan.set_defaults(func=make_manhattan)

    parser_miami = subparsers.add_parser("miami", help="Make Miami plot.")
    parser_miami.add_argument("--config-top", type=str, nargs='+', required=True, help="List of config files to plot on top.")
    parser_miami.add_argument("--config-bottom", type=str, nargs='+', required=True, help="List of config files to plot on bottom.")
    parser_miami.add_argument("--out", type=str, nargs='+', required=True, help="List of output file names.")
    parser_miami.set_defaults(func=make_miami)

    return parser.parse_args(args)


def add_coord(dfs):
    # df is assume to have "CHR" and "BP" columns. Chromosomes here are assumed to be integers.
    min_chr_bp, max_chr_bp = get_min_max_chr_bp(dfs)
    next_chrom_start = 0
    between_chr_gap = 0.005*(sum(max_chr_bp.values()) - sum(min_chr_bp.values()))
    for df in dfs:
        next_chrom_start = 0
        for chrom in sorted(min_chr_bp):
            i_chrom = df.CHR==chrom
            min_bp, max_bp = min_chr_bp[chrom], max_chr_bp[chrom]
            coord = df.loc[i_chrom, "BP"] - min_bp + next_chrom_start
            df.loc[i_chrom, "COORD"] = coord
            next_chrom_start += max_bp - min_bp + between_chr_gap
        
def get_min_max_chr_bp(dfs):
    min_max_chr_bp = []
    for df in dfs:
        df_chr_bp = df.groupby("CHR").agg({"BP":['min', 'max']}).BP
        min_max_chr_bp.append(df_chr_bp)
    df_concat = pd.concat(min_max_chr_bp, axis=1)
    min_chr_bp = df_concat.min(axis=1).to_dict()
    max_chr_bp = df_concat.max(axis=1).to_dict()
    return min_chr_bp, max_chr_bp

def drop_marginal_snps(df, p_cutoff_low, p_cutoff_high):
    # df is assume to have "P" column.
    # Add coordinates in the figure.
    # Drop non-autosomes and convert chromosomes to int.
    autosomes = [str(i) for i in range(1,23)]
    i_autosome = df.CHR.isin(autosomes)
    df = df.loc[i_autosome,:]
    i2plot = (p_cutoff_high<df.P) & (df.P<p_cutoff_low)
    df = df.loc[i2plot,:]
    print(f"{i2plot.sum()} variants will be plotted.")
    return df


def set_xlim_and_xtiks(dfs, ax):
    ax.set_xlim( (min(df.COORD.min() for df in dfs), max(df.COORD.max() for df in dfs)) )

    # add chromosome labels in the middle
    xtick_locations = []
    for chrom in range(1,23):
        bp_min = min([df.loc[df.CHR==chrom,"COORD"].min() for df in dfs if chrom in df.CHR.values])
        bp_max = min([df.loc[df.CHR==chrom,"COORD"].max() for df in dfs if chrom in df.CHR.values])
        chrom_mid = bp_min + 0.5*(bp_max - bp_min)
        xtick_locations.append(chrom_mid)

    x_tick_labels = [str(i) for i in range(1,23)]
    ax.set_xticks(xtick_locations)
    ax.set_xticklabels(x_tick_labels)


def plot_scatter(df, color1, color2, size, marker, ax, edgewidth=0, edgecolor='k', alpha=1, label=None, rasterized=True):
    color = [color1 if chrom%2==1 else color2 for chrom in df.CHR]
    _=ax.scatter(df.COORD, -np.log10(df.P), s=size, c=color, marker=marker,
                 linewidths=edgewidth, edgecolors=edgecolor, alpha=alpha, label=label, rasterized=rasterized)
    
def plot_sumstats(df, config, ax):
    plot_scatter(df.loc[df.IS_STANDARD,:], config["color1"], config["color2"], config["size"],
             config["marker"], ax, alpha=config["alpha"], label=config["legend_label"], rasterized=True)
    if config["bold"]:
        plot_scatter(df.loc[df.IS_BOLD,:], config["color1_bold"], config["color2_bold"],
                     config["size_bold"], config["marker_bold"], ax,
                     alpha=config["alpha_bold"], rasterized=True)
    if config["outlined"]:
        plot_scatter(df.loc[df.IS_OUTLINED,:], config["color1_outlined"], config["color2_outlined"],
                     config["size_outlined"], config["marker_outlined"], ax,
                     alpha=config["alpha_outlined"], edgewidth=1.5, edgecolor='k', rasterized=True)
    if config["gws_threshold"]:
        gws_thresh = -np.log10(config["gws_threshold"])
        ax.hlines([gws_thresh], 0, 1, colors='k', linestyles='dotted', transform=ax.get_yaxis_transform())


def read_sumstats(config):
    print(f"Reading sumstats from: {config['sumstats']}")
    df = pd.read_csv(config["sumstats"], delim_whitespace=True,
                     usecols=[config["chrom_col"], config["bp_col"], config["p_col"], config["id_col"]],
                     dtype={config["chrom_col"]:str, config["bp_col"]:int, config["p_col"]:float, config["id_col"]:str})
    # standardize column names
    rename_col_dict = {config["chrom_col"]:"CHR", config["bp_col"]:"BP", config["p_col"]:"P", config["id_col"]:"ID"}
    df.rename(columns=rename_col_dict, inplace=True)
    
    print(f"    {df.shape[0]} variants loaded.")
    df.rename(columns={config["id_col"]:"ID", config["chrom_col"]:"CHR", config["bp_col"]:"BP", config["p_col"]:"P"}, inplace=True)

    p_cutoff_low = float(config["p_cutoff_low"])
    p_cutoff_high = float(config["p_cutoff_high"])
    df = drop_marginal_snps(df, p_cutoff_low, p_cutoff_high)

    df["IS_OUTLINED"] = False
    if config["outlined"]:
        print(f"Reading {config['outlined']}")
        outlined = pd.read_csv(config["outlined"], header=None, names=["ID"]).squeeze()
        df["IS_OUTLINED"] = df.ID.isin(outlined)

    df["IS_BOLD"] = False
    if config["bold"]:
        print(f"Reading {config['bold']}")
        bold = pd.read_csv(config["bold"], header=None, names=["ID"]).squeeze()
        df["IS_BOLD"] = df.ID.isin(bold)
        if config["outlined"]:
            df.loc[df.IS_BOLD & df.IS_OUTLINED,"IS_BOLD"] = False

    df["IS_STANDARD"] = True
    df.loc[df.IS_BOLD | df.IS_OUTLINED, "IS_STANDARD"] = False
    return df


def make_manhattan(args):
    config_files = args.config
    configs, dfs = [], []
    for config_file in config_files:
        print(f"Reading config from: {config_file}")
        config = {}
        with open(config_file) as cf:
            exec(cf.read(),  {"__builtins__" : None}, config)
        configs.append(config)
        df = read_sumstats(config)
        dfs.append(df)
        df.CHR = df.CHR.astype(int)

    add_coord(dfs)

    fig, ax = plt.subplots(1, 1, figsize=(18,4), constrained_layout=True)    
    # configure spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_position(('outward',5))
    
    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=True,       # ticks along the bottom edge are on/off
        top=False,         # ticks along the top edge are on/off
        labelbottom=True)  # tick labels along the bottom edge are on/off
    
    ax.set_xlabel("Chromosome", fontsize=16)
    
    for df, config in zip(dfs, configs):
        plot_sumstats(df, config, ax)
        
    print("Plotting data.")
    set_xlim_and_xtiks(dfs, ax)

    ax.set_ylabel(r"$\mathrm{-log_{10}(%s)}$" % configs[-1]["y_axis_label"], fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=13)
    label_colors = [c["legend_label_color"] for c in configs]
    ax.legend(loc="upper right", fontsize=16, frameon=False, labelcolor=label_colors, markerscale=0, scatterpoints=1)

    # set y axis lims
    _, y_max = ax.get_ylim()
    y_min = min(-np.log10(df.P.max()) for df in dfs)
    ax.set_ylim((y_min, y_max))

    for outf in args.out:
        plt.savefig(outf, facecolor='w', dpi=300)
        print(f"{outf} saved.")


def make_miami(args):
    config_files = args.config_top + args.config_bottom
    configs, dfs = [], []
    for config_file in config_files:
        print(f"Reading config from: {config_file}")
        config = {}
        with open(config_file) as cf:
            exec(cf.read(),  {"__builtins__" : None}, config)
        configs.append(config)
        df = read_sumstats(config)
        dfs.append(df)
        df.CHR = df.CHR.astype(int)

    add_coord(dfs)

    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(18,8), constrained_layout=True, sharex=True)
    ax_bottom.invert_yaxis()
    ax_bottom.tick_params(axis='x', which='both', top=False, labeltop=True, bottom=False, labelbottom=False)
    ax_top.tick_params(axis='x', which='both', top=False, labeltop=False, bottom=False, labelbottom=False)
    axs = [ax_top]*len(args.config_top) + [ax_bottom]*len(args.config_bottom)
    # configure spines
    for ax in (ax_top, ax_bottom):
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_position(('outward',5))
        ax.ticklabel_format(axis='x', useOffset=False)
    
    for config, df, ax in zip(configs, dfs, axs):
        plot_sumstats(df, config, ax)
        
    set_xlim_and_xtiks(dfs, ax_bottom) # only bottom axis have tick labels on top

    ax_top.set_ylabel(r"$\mathrm{-log_{10}(%s)}$" % configs[len(args.config_top)-1]["y_axis_label"], fontsize=16)
    ax_bottom.set_ylabel(r"$\mathrm{-log_{10}(%s)}$" % configs[-1]["y_axis_label"], fontsize=16)
    ax_top.tick_params(axis='both', which='major', labelsize=13)
    ax_bottom.tick_params(axis='both', which='major', labelsize=13)
    legend_label_colors = [c["legend_label_color"] for c in configs[:len(args.config_top)]]
    ax_top.legend(loc="upper right", fontsize=16, frameon=False, labelcolor=legend_label_colors, markerscale=0, scatterpoints=1)
    legend_label_colors = [c["legend_label_color"] for c in configs[len(args.config_top):]]
    ax_bottom.legend(loc="lower right", fontsize=16, frameon=False, labelcolor=legend_label_colors, markerscale=0, scatterpoints=1)

    # set y axis lims
    _, y_max_top = ax_top.get_ylim()
    y_min_top = min(-np.log10(df.P.max()) for df in dfs[:len(args.config_top)])
    y_max_bottom, _ = ax_bottom.get_ylim()
    y_min_bottom = min(-np.log10(df.P.max()) for df in dfs[len(args.config_top):])
    if any([c["allign_y_max"] for c in configs]):
        y_max_top = max(y_max_top, y_max_bottom)
        y_max_bottom = y_max_top
    ax_top.set_ylim((y_min_top, y_max_top))
    ax_bottom.set_ylim((y_max_bottom, y_min_bottom))

    for outf in args.out:
        plt.savefig(outf, facecolor='w', dpi=300)
        print(f"{outf} saved.")


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    args.func(args)
