import pandas as pd 
from pathlib import Path 

OVERVIEW_DIR = "../SampleOverview/data"
HET_TSV = "../HetComposition/data/HetsComposition.tsv"


def import_overview(dir_path):
    """
    """
    dfs = [] 
    for fpath in Path(dir_path).glob('*xlsx'):
        # data type (raw or relative)
        type = fpath.stem.split('Classification')[0]

        df = pd.read_excel(fpath)    

        # wide to long 
        df = pd.melt(
            df, 
            id_vars=['Sample'], 
            value_vars=['Prochlorococcus', 'Synechococcus', 'Viruses', 'Other Genus', 'Unclassified'], 
            var_name='taxon_name', 
            value_name='reads' if type == 'Raw' else 'percent'
        )
        
        dfs.append(df)

    # merge raw read count with relative 
    df = pd.merge(dfs[0], dfs[1], on=['Sample', 'taxon_name'], how='outer')

    df = df.dropna(subset=['percent', 'reads'])

    return df 


def import_het(tsv_path):
    """
    """
    df = pd.read_table(tsv_path)
    df = df[['Sample', 'taxon_name', 'percent', 'reads', 'rank']]

    genus_df = df[df['rank'] == 'genus'].copy()
    class_df = df[df['rank'] == 'class'].copy()

    genus_df.drop(columns=['rank'], inplace=True)
    class_df.drop(columns=['rank'], inplace=True)

    return genus_df, class_df
   

def make_supplemental(overview_df, genus_df, class_df):
    """
    """
    with pd.ExcelWriter('data/SupplementalTable.xlsx') as writer:
        overview_df.to_excel(writer, sheet_name='Overview', index=False)
        genus_df.to_excel(writer, sheet_name='Genus', index=False)
        class_df.to_excel(writer, sheet_name='Class', index=False)

    return True 


def main():
    Path('data').mkdir(exist_ok=True)

    # process overview data 
    overview_df = import_overview(OVERVIEW_DIR)

    # process het data 
    genus_df, class_df = import_het(HET_TSV)

    # combine data into 1 supplemental table 
    make_supplemental(overview_df, genus_df, class_df)

if __name__ == "__main__":
    main()