"""
Parse bcftools query output for supplementary table. 
"""
import pandas as pd 
from pathlib import Path 
import vcf_parse_helper

samples_tsv = "../../../sample_metadata.tsv"
outdir = 'data/2.parsed_vcf'

def process_vcf(fpath, gene_df):
    """
    Returns df with mutation info from vcf file and gene info. 
    """
    df = pd.read_table(fpath, names=['CHROM', 'POS', 'REF', 'ALT', 'INDEL'])

    # determine type of mutation (SNP, insertion, deletion)
    df['mutation type'] = df.apply(
        lambda row: 'SNV' if row['INDEL'] == '.'  
        else 'deletion' if len(row['REF']) > len(row['ALT'].split(',')[0]) 
        else 'insertion',  
        axis=1
    )

    # edit Syn9451 contig name to match NCBI 
    df['CHROM'] = df['CHROM'].str.replace('contig_76', 'NZ_CP138971.1')

    # assign mutations to annotated regions
    # for rows with the same CHROM, combine every mutation row with every gene/intergenic rows 
    df = pd.merge(df, gene_df, how='cross')
    # print(gene_df['seq_id'].unique())

    # filter to keep rows where mutation is WITHIN the gene/intergenic regions and has the same chromosome 
    df = df[(df['POS'] >= df['start']) & (df['POS'] <= df['end']) & (df['seq_id'] == df['CHROM'])]

    return df

def final_formatting(df):
    """
    Edit df before saving 
    """
    df = df.rename(columns={
        'CHROM': 'contig id', 
        'POS': 'mutation position', 
        'REF': 'reference nucleotide(s)', 
        'ALT': 'alternative nucleotide(s)', 
        'start': 'gene start', 
        'end': 'gene end', 
        'ID': 'gene ID', 
        'name': 'gene name', 
        'product': 'gene product', 
    })

    df = df.drop(['INDEL', 'seq_id'], axis=1)

    df['replicate'] = df['Sample'].str.split('_').str[-1]
    df['treatment'] = df['Sample'].str.split('_').str[1]

    col_order = [
        'Sample', 
        'treatment', 
        'replicate', 
        'mutation type', 
        'mutation position', 
        'reference nucleotide(s)', 
        'alternative nucleotide(s)', 
        'gene ID', 
        'gene name', 
        'gene product', 
        'gene start', 
        'gene end', 
        'contig id', 
    ]
    df = df[col_order]

    return df 

def save_all_mutations(df, outdir):
    """
    """
    # group by experiment (Syn / M4A1A / M4XE)
    df['exp'] = df['Sample'].str.split('_').str[0]
    groups = df.groupby(['exp'], sort=False)

    output_xlsx_path = f"{outdir}/All-Mutations.xlsx"
    with pd.ExcelWriter(output_xlsx_path) as file:
        for index, gdf in groups:
            exp = index[0]

            gdf = gdf.drop(['exp'], axis=1)

            gdf.to_excel(file, sheet_name=exp, index=False)

    print(f'DataFrame of all mutations saved to {output_xlsx_path}!')



def filter_mutations(df):
    """
    Returns df with mutation filtering: 
        1. Condense mutations by genic/intergenic 
        2. Remove mutations present in < 2 bio rep 
        2. Obtain mutations present in: 
            a. NPC only 
            b. both 1st and 2nd infection 
    """
    # condense mutation in same genic/intergenic region into 1 row 
    # (e.g. '1 insertion, 3 SNV' if region has both)
    df['mutation type'] = (
        df.groupby(['treatment', 'replicate', 'gene ID'])['mutation type']
        .transform(lambda x: ', '.join(f"{count} {value}" for value, count in x.value_counts().items()))
    )
    df = df.drop_duplicates(subset=['Sample', 'gene ID'], keep='first')

    # REP filtering: group by treatment and gene/intergenic region --> remove mutations in < 2 reps 
    df = df.groupby(['treatment', 'gene ID']).filter(lambda x: x['replicate'].nunique() > 1)

    # TREATMENT filtering: group by gene/intergenic region
    # keep mutations in: NPC only or 1st AND 2nd infection 
    infection_treatments = {
        x for x in df['treatment'].unique() if x!='NPC'
    }

    df = df.groupby(['gene ID']).filter(
        lambda x: (
            set(x['treatment']) == {'NPC'} or 
            set(x['treatment']) == infection_treatments
        )
    )    

    return df

def supplemental_filter(df, outdir):
    """
    """
    # group by experiment (Syn / M4A1A / M4XE)
    df['exp'] = df['Sample'].str.split('_').str[0]
    groups = df.groupby(['exp'], sort=False)

    output_xlsx_path = f"{outdir}/Filtered-Mutations.xlsx"
    with pd.ExcelWriter(output_xlsx_path) as file:
        for index, gdf in groups:
            exp = index[0]

            gdf = gdf.drop(['exp'], axis=1)

            gdf = filter_mutations(gdf)

            gdf.to_excel(file, sheet_name=exp, index=False)

    print(f'DataFrame of all mutations saved to {output_xlsx_path}!')


def plot_format(df, outdir):
    """
    Saves df for plotting. 
    """
    # group by experiment and gene id 
    groups = df.groupby(['Sample', 'gene ID'], sort=False)

    data = []
    for index, gdf in groups:
        group_data = {
            'Sample': gdf['Sample'].iloc[0], 
            'gene ID': gdf['gene ID'].iloc[0],
            'mutation count': len(gdf), 
            'gene midpoint': int((gdf['gene end'].iloc[0] + gdf['gene start'].iloc[0]) / 2), 
            'contig id': gdf['contig id'].iloc[0],
        }

        data.append(group_data)

    df = pd.DataFrame(data)

    # save plotting data 
    df.to_csv(f'{outdir}/plot_table.tsv', sep='\t', index=False)

def main():
    # make output dir 
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # import gff files 
    ref_genome_dir = "../../ref_genomes"
    gene_df = vcf_parse_helper.obtain_genes(ref_genome_dir)

    # import samples tsv
    sdf = pd.read_table(samples_tsv)

    dfs = []
    for fpath in Path('data/1.bcftools_out').glob('*'):
        fname = fpath.stem.split('_ref')[0]  # obtain sample name 
        df = process_vcf(fpath, gene_df)
        df['FileName'] = fname

        dfs.append(df)
        
    df = pd.concat(dfs)
    df = pd.merge(sdf, df, on=['FileName'], how='inner')

    # format df before saving 
    df = final_formatting(df)  

    # save all mutations 
    save_all_mutations(df, outdir)  

    # filter for supplemental table 
    supplemental_filter(df, outdir)  

    # format for plotting
    plot_format(df, outdir)   


main()