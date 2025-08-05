"""
"""
import pandas as pd
import numpy as np
from pathlib import Path 
import matplotlib.pyplot as plt

SAMPLES_TSV = "../../../sample_metadata.tsv"

def assign_taxon_colors(df):
    """
    """
    # group by rank 
    groups = df.groupby(['rank'], sort=False)

    color_dict_list = []

    for index, gdf in groups:
        rank = index[0]
        if rank == 'class':
            cmap = plt.get_cmap('tab20')
        else:
            cmap = plt.get_cmap('tab20')

        # obtain unique taxon
        unique_taxon = gdf['taxon_name'].unique()

        # map each unique value to a color 
        colors = [cmap(i / len(unique_taxon)) for i in range(len(unique_taxon))]
        taxon_to_color = dict(zip(unique_taxon, colors))

        color_dict_list.append(taxon_to_color)

    color_dict = color_dict_list[0] | color_dict_list[1]

    # save file 
    with open('data/HetsColors.tsv', 'w') as file:
        file.write('taxon_id\tcolor\n')
        for key, value in color_dict.items():
            value = tuple([float(x) for x in value])
            file.write(f'{key}\t{value}\n')

    return color_dict

def plot(df, sdf, plot_name, color_dict):
    """
    """
    # reformat (pivot) df 
    df = df.pivot(index='Sample', columns='taxon_name', values='percent').reset_index()

    # reorder samples 
    df["Sample"] = pd.Categorical(df["Sample"], categories=sdf["Sample"], ordered=True)
    df = df.sort_values("Sample").reset_index(drop=True)

    # obtain bar order
    order_list = df['Sample'].values.tolist()
    order_list.reverse()
    df = df.set_index('Sample').reindex(order_list)

    # obtain colors 
    colors = [color_dict[taxon] for taxon in df.columns.tolist()]

    # plotting 
    fig, ax = plt.subplots(1, 1, figsize=(12, 4)) # (w, h)
    df.plot.barh(stacked=True, ax=ax, width=0.8, edgecolor='black', linewidth=0.7, color=colors)

    # move legend to the right 
    ax.legend(
        loc='center left', bbox_to_anchor=(1, 0.5), 
        fontsize=13, frameon=False,  # remove square frame
        title=f'', title_fontsize=16, 
        handleheight=1, handlelength=1,  # make markers square instead of rectangles
        handletextpad=0.5, labelspacing=0,  # no spacing between markers vertically 
        edgecolor='0',
        prop={'style': 'italic'}, # italicize legend labels 
    )

    # leave space for legend and y-axis tick labels 
    plt.subplots_adjust(left=0.3, right=0.6)

    # titles
    ax.set_ylabel('')
    # ax.set_title(plot_name)
    plt.suptitle(plot_name)

    # plt.tight_layout()
    plt.savefig(f'data/{plot_name}.png')
    plt.savefig(f'data/{plot_name}.svg')
    plt.close()

def main():
    Path('data').mkdir(exist_ok=True, parents=True)

    df = pd.read_table('data/HetsComposition.tsv')
    sdf = pd.read_table(SAMPLES_TSV)

    # obtain colors 
    color_dict = assign_taxon_colors(df)

    # group for plots 
    plot_groups = df.groupby(['Treatment', 'rank'], sort=False)

    for index, df in plot_groups:
        plot_name = f'{index[0]}_{index[1]}'
        plot(df, sdf, plot_name, color_dict)


main()