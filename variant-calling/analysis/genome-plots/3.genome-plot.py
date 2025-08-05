"""
"""
import pandas as pd 
from pathlib import Path 
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

def obtain_genome_len(dir_path):
    """
    Returns df with length of MED4 & WH8102 reference genomes. 
    """
    data = {'seq_id': [], 'genome_len': []}
    for genome_fpath in Path(dir_path).glob('*fasta'):
        # obtain seq_id (contig name)
        with open(genome_fpath, 'r') as file:
            seq_id = file.readline().strip().split(' ')[0].replace('>', '')

        # obtain genome length 
        genome_len = 0

        for record in SeqIO.parse(genome_fpath, "fasta"):
            genome_len += len(record.seq)

        data['seq_id'].append(seq_id)
        data['genome_len'].append(genome_len)

    df = pd.DataFrame(data)

    return df

def format_tick_labels(x, pos):
    """
    Format tick labels to display thousands as 'K' and millions as 'M'
    """
    if x >= 1000000:
        # Format as millions
        return f'{x / 1000000:.1f}M' if x % 1000000 != 0 else f'{int(x / 1000000)}M'
    elif x >= 1000:
        # Format as thousands
        return f'{int(x / 1000)}K'
    return str(int(x))

def plot_mut_on_genome(df, exp, genome_len_df):
    """
    """
    # obtain genome length from contig name (seq_id)
    contig_id = df[df['exp'] == exp]['contig id'].iloc[0]
    genome_len = genome_len_df[genome_len_df['seq_id'] == contig_id]['genome_len'].iloc[0]

    # obtain color intensity (based on number of mutation at each gene)
    color_norm = mcolors.Normalize(vmin=df['mutation count'].min(), vmax=df['mutation count'].max())

    # split dfs into reps
    groups = df.groupby(['treatment', 'rep'], sort=False)

    # plotting 
    fig, axes = plt.subplots(len(groups), 1, figsize=(15, 5), sharex=True)
    axes = axes.flatten() 

    for ax_index, ((treatment, rep), ax_df) in enumerate(groups):
        # obtain axes to plot data on 
        ax = axes[ax_index]
        ax_title = f'{treatment} {rep}'

        # 1. plot horizontal bar (len of genome)
        ax.barh([0], genome_len, color='lightgray', alpha=0.1)  

        # 2. plot each mutated gene 
        for gene_midpoint, region_type, mut_count in zip(ax_df['gene midpoint'], ax_df['region type'], ax_df['mutation count']):
            base_color = 'red' if region_type == 'Gene' else 'blue'

            # Adjust alpha by mutation count
            adjusted_color = mcolors.to_rgba(base_color, alpha=color_norm(mut_count))

            # plot line of mutated region 
            ax.axvline(x=gene_midpoint, color=base_color, linestyle='-', linewidth=0.5)


        # other axes edits 
        ax.set_xlim(-2, genome_len+2)  # Set x-axis limit from -1 to genome_len+1
        ax.set_ylabel(ax_title, rotation=0, ha='right', va='center')   # va='top', , labelpad=50

        if ax_index != (len(axes) - 1):
            # Hide the ticks for all but last axes' axis
            ax.tick_params(axis='x', direction='in', length=0) 
        else:
            # edit x-ticks on the last axes' axis  
            tick_interval = 100000
            ticks = range(0, genome_len + 1, tick_interval)
            ax.set_xticks(ticks)
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_tick_labels))

        # remove y-axis-ticks for all axes
        ax.set_yticks(()) 
        ax.set_yticklabels([])  

    # add legend for colors
    legend_elements = [
        Line2D([0], [0], color='red', lw=2, label='Gene'),
        Line2D([0], [0], color='blue', lw=2, label='Intergenic')
    ]

    fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))

    plt.suptitle(f'{exp} Mutated Regions')
    plt.tight_layout()

    outdir = "data/3.genome_plots"
    Path(outdir).mkdir(parents=True, exist_ok=True)
    plt.savefig(f'{outdir}/{exp}.png', bbox_inches='tight')
    plt.savefig(f'{outdir}/{exp}.svg', format='svg', bbox_inches='tight')



def plot_vcf(fpath, genome_len_df):
    """
    """
    df = pd.read_table(fpath)
    
    df['exp'] = df['Sample'].str.split('_').str[0]
    
    df['treatment'] = df['Sample'].str.split('_').str[1]
    df['rep'] = df['Sample'].str.split('_').str[-1]
    df['region type'] = df['gene ID'].apply(lambda x: 'Intergenic' if 'intergenic' not in x else 'Gene')

    # group by experiment and plot each 
    for index, gdf in df.groupby(['exp']):
        plot_mut_on_genome(gdf, index[0], genome_len_df)

def main():
    # obtain genome lengths (MED4 & S9451)
    genome_dir = "../../ref_genomes"
    genome_len_df = obtain_genome_len(genome_dir)

    # mutation plotting 
    parsed_df = "data/2.parsed_vcf/plot_table.tsv"
    plot_vcf(parsed_df, genome_len_df)

main()