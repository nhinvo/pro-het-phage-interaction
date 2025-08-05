"""
Purpose: plot all positions in cmt table. 
    - Note: these positions are from variant calling, have not been filtered. 
"""
import yaml
import pandas as pd 
from pathlib import Path 
import matplotlib.pyplot as plt

SAMPLES_TSV = "../../../sample_metadata.tsv"

def load_yaml(yaml_path):
    """
    Returns contents of a yaml file as a dictionary. 
    """
    with open(yaml_path, 'r') as file:
        yaml_content = yaml.safe_load(file)

    return yaml_content

def extract_lowercase(string):
    """
    Returns lowercase letter in row value only. 
        - e.g. 'parent' from 'parentL', 'adh' from 'adh2'. 
    """
    return ''.join([char for char in string if char.islower()])

def extract_number(string):
    """
    Returns numbers in row value only. 
        - e.g. '2' from 'adh2', '12' from 'liq12'
        - For parent samples: parentL = '0', parentA = '1'. 
    """
    # obtain numerical values in sample string 
    number = ''.join([char for char in string if char.isdigit()])

    # parent samples have no numbers = empty string
    if number == '':
        # obtain the parent type ('L' or 'A')
        parent_type = ''.join([char for char in string if char.isupper()])
        # assign number based on type 
        number = '0' if 'L' in parent_type else '1'

    return number 

def import_cmt_df(fpath, sdf):
    """
    Import and process cmt table. 
    """
    df = pd.read_table(fpath)
    df = df.rename(columns={'sample': 'FileName'})
    df = pd.merge(df, sdf, on=['FileName'], how='left')

    # obtain treatment & rep from sample name
    df['Treatment'] = df['Sample'].str.split('_').str[1]
    df['Replicate'] = df['Sample'].str[-1]

    # save table
    # group = fpath.stem.split('_')[1]
    # Path('data').mkdir(exist_ok=True, parents=True)
    # df.to_csv(f'data/sorted_pos_data.tsv', sep='\t', index=False)

    return df 

def plot_bar_position(df, position, dataset):
    """
    For position df (base count for 1 position), plot count of each base 
    in the forward and reverse direction, then save plot. 
    """
    # split into 2 df - fwd and rev count 
    fwd_counts = df[['Sample', 'A', 'C', 'G', 'T']].set_index('Sample').fillna(0)
    rev_counts = df[['Sample', 'a', 'c', 'g', 't']].set_index('Sample').fillna(0)

    # plotting 
    fig, ax = plt.subplots(figsize=(7, 5))
    fwd_counts.plot(kind='bar', stacked=True, width=0.4, ax=ax, position=1)
    rev_counts.plot(kind='bar', stacked=True, width=0.4, ax=ax, position=0, legend=None)

    # axes labels
    ax.set_title(f'Base count per sample for position {position}')
    ax.set_ylabel('Read Count')
    ax.set_xlabel('Sample')

    # bold reference on legend 
    ref_base = str(df['ref'].iloc[0])
    legend_dict={'A':'A', 'C':'C', 'G':'G', 'T':'T'}
    legend_list=[]  # list of legend texts to be plotted 
    legend_dict[ref_base]=r'$\bf{{letter}}$'.replace('letter', ref_base)  # latex formatting for bolded reference 
    for key, value in legend_dict.items():
        legend_list.append(value)
    ax.legend(legend_list, loc='upper left')

    plt.tight_layout()

    # save plot
    png_outdir = f'data/position_barplots/{dataset}'
    # svg_outdir = f'data/1.position_barplots/svg'

    Path(png_outdir).mkdir(parents=True, exist_ok=True)
    # Path(svg_outdir).mkdir(parents=True, exist_ok=True)

    plt.savefig(f'{png_outdir}/{position}.png')
    # plt.savefig(f'{svg_outdir}/{position}.svg')

    plt.close()   # clear plot 


def plot_bars(df, dataset):
    """
    For each position in each dataset, plot base calls for all samples. 
    """
    # group by position and plot base call at that pos for all samples
    pos_groups = df.groupby(['position'])
    for index, pos_df in pos_groups:
        position = index[0]
        plot_bar_position(pos_df, position, dataset)

def main():
    sdf = pd.read_table(SAMPLES_TSV)

    cmt_dir = "data/parsed_cmt"
    for fpath in Path(cmt_dir).glob('*tsv'):
        dataset = fpath.stem.split('_')[1]

        # prep cmt file 
        cmt_df = import_cmt_df(fpath, sdf)

        # plot all positions in cmt 
        plot_bars(cmt_df, dataset)
    
if __name__ == "__main__":
    main()