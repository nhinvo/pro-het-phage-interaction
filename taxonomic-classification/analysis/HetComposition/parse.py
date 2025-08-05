import os 
import pandas as pd 
from pathlib import Path 

MIN_PCNT = 1
SAMPLES_TSV = "../../../sample_metadata.tsv"
KAIJU_SUMMARY_DIR = "../../modified-prosyntax-workflow/scratch/kaiju_summary"

def parse_summary_file(fpath):
    """
    Import kaiju summary output file. 
    """
    # metadata from file path 
    rank = fpath.stem.split('_')[-1]
    sample = '_'.join(fpath.stem.split('_')[0:2])

    # import & edit table 
    df = pd.read_table(fpath)
    df['rank'] = rank
    df['FileName'] = sample
    df = df.drop(columns=['file'])
    
    return df 

def process_summary_data(kaiju_summary_dir, sdf):
    """
    """
    out_fpath = 'data/parsed_summary.tsv'

    # check if the output file already exists
    if os.path.exists(out_fpath):
        df = pd.read_table(out_fpath)
        return df

    # import and obtain table if not 
    dfs = []
    for fpath in Path(kaiju_summary_dir).glob('*tsv'):
        df = parse_summary_file(fpath)
        dfs.append(df)

    df = pd.concat(dfs)
    df = pd.merge(sdf, df, on=['FileName'], how='inner')
    df.to_csv(out_fpath, sep='\t', index=False)

    return df

def calculate_pcnt(df, cols):
    """
    Calculate percentage (reads) for each group.
    """
    groups = df.groupby(cols, sort=False)

    dfs = []
    for index, gdf in groups:
        total_read = gdf['reads'].sum()
        gdf['percent'] = (gdf['reads'] / total_read) * 100 
        dfs.append(gdf)

    df = pd.concat(dfs)

    return df 

def pcnt_grouping(df, pcnt):
    """
    """
    # group before processing 
    groups = df.groupby(['Treatment', 'rank'], sort=False)

    dfs = []

    for index, gdf in groups:
        # obtain taxon rows that passed threshold
        taxon_passed = gdf[gdf['percent'] > float(pcnt)]['taxon_name'].values.tolist()
        taxon_passed = set(taxon_passed)

        # group by sample and combine taxon that did not pass threshold
        sample_groups = gdf.groupby(['FileName'], sort=False)

        for index, sdf in sample_groups:
            # filter for taxon that passed  
            genus_sdf = sdf[sdf['taxon_name'].isin(taxon_passed)]

            # combine remaining taxons 
            other_sdf = sdf[~sdf['taxon_name'].isin(taxon_passed)]
            new_row = pd.DataFrame({
                'taxon_name': ['Others'], 
                'reads': [other_sdf['reads'].sum()],  
                'percent': [other_sdf['percent'].sum()], 
                'FileName': [index[0]], 
                'Sample': other_sdf['Sample'].iloc[0],
                'Treatment': other_sdf['Treatment'].iloc[0],
                'rank': other_sdf['rank'].iloc[0],
            })

            # create new df 
            sdf = pd.concat([
                genus_sdf, 
                new_row
            ])

            dfs.append(sdf)

    df = pd.concat(dfs)

    return df

def excel_save(df):
    """
    Save df as excel file. 
    """
    # group data into tabs 
    tabs = df.groupby(['Treatment', 'rank'], sort=False)

    with pd.ExcelWriter('data/HetsComposition.xlsx') as writer:
        for index, df in tabs:
            tab_name = f'{index[0]}_{index[1]}'
            
            col_order = ['Sample', 'taxon_name', 'reads', 'percent']
            df = df[col_order]
            
            df.to_excel(writer, sheet_name=tab_name, index=False)

def process_data(df):
    """
    """
    # filter for data of interest
    df = df[df['rank'].str.contains('genus|class')].copy()

    # remove certain taxons
    remove = ['Cyanophyceae', 
              'Prochlorococcaceae',
              'Synechococcaceae',
              'Prochlorococcus',
              'Synechococcus',
              'Synechococcales',
              'Cyanobacteriota',
              'unclassified', 
              'cannot be assigned to a',
              'Viruses',
              'Thermus', 
              ]
    remove = '|'.join(remove)
    df = df[~df['taxon_name'].str.contains(remove)].copy()

    # re-calculate percentages
    cols = ['Sample', 'rank']
    df = calculate_pcnt(df, cols)

    # for each dataset, obtain a list of taxon that is above 1% in at least 1 sample
    # group the remaining into "Others"
    df = pcnt_grouping(df, MIN_PCNT)

    # save tsv data 
    df.to_csv('data/HetsComposition.tsv', sep='\t', index=False)
    excel_save(df)


def main():
    # create output dir 
    Path('data').mkdir(exist_ok=True, parents=True)  

    # import samples tsv
    sdf = pd.read_table(SAMPLES_TSV)
    sdf['Treatment'] = sdf['Sample'].str.split('_').str[0]

    # import summary df
    df = process_summary_data(KAIJU_SUMMARY_DIR, sdf)

    # process summary data for plotting 
    process_data(df)


if __name__ == "__main__":
    main()