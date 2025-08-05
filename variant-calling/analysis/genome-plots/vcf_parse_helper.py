"""
Purpose: script for parsing gff file for merging with gff. 
"""
import pandas as pd 
from pathlib import Path 

def create_attr_dict(attr_row):
    """
    Function from Konnor. 
    """
    d = {}
    for key_value_pair in attr_row.split(";"):
        k_v_list = key_value_pair.split(sep="=", maxsplit=1)
        if len(k_v_list) == 2:
            k, v = k_v_list
            d[k] = v
    return d

def generate_gff_df(gff_file):
    """
    Function from Konnor. 
    """
    ## This method is heavily influenced from GFFpandas
    df = pd.read_csv(gff_file, sep="\t", comment="#",
            names=[
                "seq_id",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ],
    )

    attr_dict_series = df.attributes.apply(
        lambda attributes: create_attr_dict(attributes)
    )

    key_set = set()
    attr_dict_series.apply(lambda at_dic: key_set.update(list(at_dic.keys())))
    for attr in sorted(list(key_set)):
        df[attr] = attr_dict_series.apply(
            lambda attr_dict: attr_dict.get(attr)
        )

    return df

def add_intergenic(df):
    """
    Create rows for intergenic regions. 
    Note: if gff file has more than 1 contig ('>'), function will be wrong
        as it uses the first seq_id for new intergenic rows. 
    """
    seq_id = df['seq_id'].unique()[0]
    intergenic_rows = []

    for i in range(len(df) - 1):
        # obtain rows 
        current_row = df.iloc[i]
        next_row = df.iloc[i + 1]

        # check if there is gap (intergenic) regions between rows 
        if next_row['start'] > current_row['end'] + 1:
            # create new row for this intergenic region 
            new_row = {
                'seq_id': seq_id, 
                'start': current_row['end'] + 1,  # Start right after the 'end' of the current row
                'end': next_row['start'] - 1,  # End right before the 'start' of the next row
                'ID': f'intergenic-{current_row['end'] + 1}:{next_row['start'] - 1}',
                'name': '',
                'product': ''
            }
            # Append the new row to the list
            intergenic_rows.append(new_row)

    intergenic_df = pd.DataFrame(intergenic_rows)

    # add intergenic to original df 
    df = pd.concat([df, intergenic_df]).sort_values(by='start').reset_index(drop=True)

    return df

def obtain_genes(dir_path):
    """
    Returns df of both MED4 and WH8102 with:
        - chromosome id
        - start and end coords for genes (CDS) and intergenic regions
        - gene annotation 
    """
    dfs = []  
    for gff_fpath in Path(dir_path).glob('*.gff'):
        genome_name = gff_fpath.name.split('_')[0]  # i.e. MED4 / MITS9451 

        # import gff 
        df = generate_gff_df(gff_fpath)

        # filter rows 
        if 'MITS9451' in genome_name:
            df = df[df['type'] == 'CDS'].copy()
            gene_col = 'Name'

        else:
            gene_col = 'em_Preferred_name'

        # df edits
        df = df[['seq_id', 'start', 'end', 'ID', gene_col, 'product']]  
        df = df.rename(columns={
            gene_col: 'name', 
        })

        # create intergenic rows 
        df = add_intergenic(df)

        # add to list of gff dfs 
        dfs.append(df)
    
    # combine MED4 and WH8102 gff dfs
    df = pd.concat(dfs)

    return df