"""
Purpose: to copy read files into this input dir. 
"""
import shutil 
import pandas as pd
from pathlib import Path 

def copy_files():
    """
    Copy raw read files into raw_reads/ directory.
    """
    # import df of samples and their paths 
    df = pd.read_table('data/sample_paths.tsv')

    for index, row in df.iterrows():
        dir_name = row['sample']

        print(f"Copying files for {dir_name}...")

        # make dir to copy input over if does not exist 
        output_dir = f'raw_reads/{dir_name}'
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # obtain read files and copy over 
        fwd_read = row['forward read']
        rev_read = row['reverse read']

        for fpath in [fwd_read, rev_read]:
            try: 
                shutil.copy(fpath, f'{output_dir}/{fpath.split("/")[-1]}')
            except Exception as e:
                print(f'An error occured: {e}')

    return True 

def make_samples_csv():
    """
    """
    raw_read_realpath = "/orcd/data/chisholm/002/hstor001/nvo/projects/25_AllisonPhageVariantCalling/variant-calling/modified-widevariant-workflow/0-input/raw_reads"
    df = pd.read_table('data/sample_metadata.tsv')

    # wide variant cols 
    df['Group'] = df['Sample'].str.split(' ').str[0]
    df['Reference'] = df['Group'].str.split('-').str[0]
    df['Reference'] = df['Reference'].str.replace('Xe', '')
    df['Outgroup'] = 0
    df['Sample'] = df['FileName']
    df['Path'] = f"{raw_read_realpath}/" + df['Sample'] + '/'

    # reorder cols
    df = df[['Path', 'Sample', 'FileName', 'Reference', 'Group', 'Outgroup']]
    df.to_csv('data/samples.csv', sep=',', index=False)

    return True 

def main():
    # copy raw read files 
    copy_files()

    # make samples.csv for widevariant input 
    make_samples_csv()


if __name__ == "__main__":
    main()