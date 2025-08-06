import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path 

def obtain_bar_order(df):
    """
    Bar order. 
    """
    # obtain sample names and reverse order for bar plot 
    order_list = df['Sample'].values.tolist()
    order_list.reverse()

    return order_list

def plot(df, treatment, data_type):
    """
    """
    # col filtering 
    df.drop('Treatment', axis=1, inplace=True)
    if data_type == 'Raw':
        # remove "Total" from raw read count 
        df.drop('Total', axis=1, inplace=True)

    if 'S9451' in treatment:
        # remove pro col from syn samples 
        df.drop('Prochlorococcus', axis=1, inplace=True)
    else:
        # remove syn col from pro samples 
        df.drop('Synechococcus', axis=1, inplace=True)

    # plotting 
    fig, ax = plt.subplots(1, 1, figsize=(11, 4))

    # order for samples on plot 
    bar_order = obtain_bar_order(df)
    df = df.set_index('Sample')
    df = df.reindex(bar_order)

    # plot bar
    taxon_color = '#fa9fb5' if 'S9451' in str(treatment) else '#74c476'
    colors = [taxon_color, '#fe9929', '#6baed6', '#b5b5b5']
    df.plot.barh(stacked=True, ax=ax, width=0.8, edgecolor='black', linewidth=0.7, color=colors)

    ax.legend(
        loc='center left', bbox_to_anchor=(1, 0.5), 
        fontsize=13, frameon=False,  # remove square frame
        title=f'', title_fontsize=16, 
        handleheight=1, handlelength=1,  # make markers square instead of rectangles
        handletextpad=0.5, labelspacing=0,  # no spacing between markers vertically 
    )

    # leave space for legend and y-axis tick labels 
    plt.subplots_adjust(left=0.3, right=0.6)

    if data_type == 'Raw':
        # put raw read count on log scale y-axis 
        plt.xscale('log')

    ax.set_ylabel('')
    ax.set_xlabel('')
    plt.suptitle(f'{treatment} Sample Composition')

    # plt.tight_layout()
    plt.savefig(f'data/{treatment}_Composition_{data_type}.png')
    plt.savefig(f'data/{treatment}_Composition_{data_type}.svg')
    plt.close()

    

def main():
    for fpath in Path('data').glob('*Summary.xlsx'):
        # "Raw" or "Relative"
        data_type = fpath.stem.replace('ClassificationSummary', '')

        # import df and group by treatment 
        df = pd.read_excel(fpath)
        groups = df.groupby(['Treatment'], sort=False)

        # plot each treatment 
        for index, df in groups:
            treatment = index[0]
            plot(df, treatment, data_type)


main()