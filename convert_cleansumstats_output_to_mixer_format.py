import pandas as pd
import numpy as np
import sys

# python convert_cleansumstats_output_to_mixer_format.py /cluster/projects/p697/projects/SUMSTATv3/v3.1/STD_GRCh37/CTG_COG_2018.sumstats.gz CTG_COG_2018_mixer.sumstats.gz
if __name__ == "__main__":
    fname = sys.argv[1]
    fname_out = sys.argv[2]
    print(f'processing {fname} -> {fname_out}...')
    df = pd.read_csv(fname, sep='\t', dtype=str)
    idx = df['CHR'].astype('str').str.lower().str.replace('chr', '').isin([str(i) for i in range(1, 23)])
    print(f'keep autosomes only (CHR column contains 1-22): {np.sum(~idx)} variants removed, {np.sum(idx)} variants retained')
    df = df[idx].copy()
    df['CHR'] = df['CHR'].astype(int)
    print('original columns: ' + ' '.join(df.columns))
    if 'POS' in df.columns: df.rename(columns={'POS':'BP'}, inplace=True)
    if 'RSID' in df.columns: df.rename(columns={'RSID':'SNP'}, inplace=True)
    if 'EffectAllele' in df.columns: df.rename(columns={'EffectAllele':'A1'}, inplace=True)
    if 'OtherAllele' in df.columns: df.rename(columns={'OtherAllele':'A2'}, inplace=True)
    if 'B' in df.columns: df.rename(columns={'B':'BETA'}, inplace=True)
    if 'EAF' in df.columns: df.rename(columns={'EAF':'FRQ'}, inplace=True)
    print('renamed  columns: ' + ' '.join(df.columns))
    print(f'writing {fname_out}...')    
    df.to_csv(fname_out, sep='\t', index=False)
    print('done.')
