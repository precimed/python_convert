# Download reference data from http://ctg.cncr.nl/software/magma (for example g1000_eur)
# Then you can run the tool as follows:
#    python make_maf_vector.py --ref 2558411_ref.bim --bfile merged --savemat mafvec.mat


import pandas as pd
import numpy as np
import argparse
import sys
import scipy.io as sio


def parse_args(args):
    parser = argparse.ArgumentParser(description="Generate LD matrix from genotype matrix")
    parser.add_argument("--ref", type=str, help="Reference file (for example 2558411_ref.bim or 9279485_ref.bim.")
    parser.add_argument("--bfile", type=str, help="Genotypes in plink binary format (only bim and frq files are required)")
    parser.add_argument("--savemat", default=None, type=str, help="Generate matfile for Matlab.")
    return parser.parse_args(args)

def make_maf_vector(ref, bfile, savemat):
    print('Reading {0}...'.format(ref))
    ref = pd.read_csv(ref, delim_whitespace=True)
    nsnp = ref.shape[0]
    chrpos_to_id = dict([((chr, pos), index) for chr, pos, index in zip(ref['CHR'], ref['BP'], ref.index)])
    if len(chrpos_to_id) != nsnp: raise ValueError("Duplicated CHR:POS pairs found in the reference file")

    print('Reading {0}...'.format(bfile + '.bim'))
    df_bim = pd.read_csv(bfile + '.bim', delim_whitespace=True, header=None)
    df_bim.columns=['CHR','SNP','GP','POS','A1','A2']
    print('Reading {0}...'.format(bfile + '.frq'))
    df_frq = pd.read_csv(bfile + '.frq', delim_whitespace=True)
    df_frq['POS'] = df_bim['POS']  # Assume that bim and frq files are aligned
    df_frq['INDEX'] = [chrpos_to_id.get((chr, pos), -1) for chr, pos in zip(df_frq['CHR'], df_frq['POS'])]

    mafvec = np.zeros((nsnp, 1)); mafvec[:] = np.NAN
    df = df_frq[df_frq['INDEX'] != -1]
    for index, value in zip(df['INDEX'], df['MAF']):
        mafvec[index] = value

    print('Found {0} SNPs with non-zero MAF, {1} with zero MAF, {2} missing in genotypes --- {3} SNPs in total'.format(
          (~np.isnan(mafvec) & (mafvec>0)).sum(), (mafvec == 0).sum(), np.isnan(mafvec).sum(), len(mafvec)))

    print('Saving result to {0}...'.format(savemat))
    sio.savemat(savemat, {'mafvec':mafvec}, format='5', do_compression=False, oned_as='column')


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    make_maf_vector(ref=args.ref, bfile=args.bfile, savemat=args.savemat)
    print("Done.")
