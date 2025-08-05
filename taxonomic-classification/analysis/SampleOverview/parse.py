"""
"""
import pandas as pd 
from pathlib import Path 

SUMMARY_FPATH = "../../modified-prosyntax-workflow/results/summary_read_count.tsv"
SAMPLES_TSV = "../../../sample_metadata.tsv"

def import_raw_count(fpath, sdf):
    """
    Import raw count, save df, and return raw count df for further parsing 
    """
    # import raw count df 
    df = pd.read_table(fpath)
    df = df[df['summary_type'] != 'percent'].copy()
    df = df.rename(columns={'sample_name':'FileName'})

    # merge with samples df 
    df = pd.merge(df, sdf, on=['FileName'], how='outer')

    # obtain ref genome 
    df['Group'] = df['Sample'].str.split(' ').str[0]

    groups = df.groupby(['Group'])

    dfs = []
    long_dfs = []
    for index, df in groups:
        group = index[0]

        # obtain taxon to remove or keep from final table 
        remove_taxon = 'Prochlorococcus' if 'MITS9451Xe' in str(group) else 'Synechococcus'
        keep_taxon = 'Synechococcus' if 'MITS9451Xe' in str(group) else 'Prochlorococcus'

        # classification str edit 
        df = df.rename(columns={
            'other_genus': 'Other Genus', 'unclassified': 'Unclassified', 
        })

        # add pro/syn to unclassified (depending on dataset)
        df['Unclassified'] = df['Unclassified'] + df[remove_taxon]
        
        # col filtering - remove Pro/Syn counts from list
        df = df[['Sample', keep_taxon, 'Viruses', 'Other Genus', 'Unclassified']]

        # undo the pivot for relative calculation later 
        long_df = pd.melt(df, id_vars=['Sample'], var_name='taxon_name', value_name='reads')
        long_dfs.append(long_df)

        # obtain total counts 
        df['Total'] = df[keep_taxon] + df['Viruses'] + df['Other Genus'] + df['Unclassified']

        dfs.append(df)

    # save raw data
    df = pd.concat(dfs)
    df['Treatment'] = df['Sample'].str.split(' ').str[0]
    df = df[['Sample', 'Treatment', 'Prochlorococcus', 'Synechococcus', 'Viruses', 'Other Genus', 'Unclassified', 'Total']]
    df.to_excel(f'data/RawClassificationSummary.xlsx', index=False)

    # concat long_df for relative count 
    long_df = pd.concat(long_dfs)

    # return the raw read count long format for futher parsing 
    return long_df

def calculate_relative(df):
    """
    Obtain relative abundance. 
    """
    # group by sample 
    sample_groups = df.groupby('Sample')

    dfs = []
    for index, sdf in sample_groups:
        total_count = sdf['reads'].sum()
        sdf['percent'] = (sdf['reads'] / total_count) * 100
        dfs.append(sdf)

    df = pd.concat(dfs)
    df = df.pivot(index=['Sample'], columns='taxon_name', values='percent').reset_index()   
    df['Treatment'] = df['Sample'].str.split(' ').str[0]
    df = df[['Sample', 'Treatment', 'Prochlorococcus', 'Synechococcus', 'Viruses', 'Other Genus', 'Unclassified']]
    
    df.to_excel(f'data/RelativeClassificationSummary.xlsx', index=False)

    return True 

def process_summary(summary_fpath, samples_tsv):
    """
    """
    # samples df 
    sdf = pd.read_table(samples_tsv)

    # import raw count df 
    df = import_raw_count(summary_fpath, sdf)

    # obtain relative count df 
    calculate_relative(df)



def main():
    Path('data').mkdir(exist_ok=True)

    process_summary(SUMMARY_FPATH, SAMPLES_TSV)



main()