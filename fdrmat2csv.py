### CONVERT PLEIOFDR RESULTS TO CSV ########################

# -- Modules -------------------------

import pandas as pd
import scipy.io as sio
import os
import argparse

# -- Parse arguments -------------------------

parser = argparse.ArgumentParser(description="Convert PleioFDR result.mat file to csv")
requiredNamed = parser.add_argument_group('Required arguments')
requiredNamed.add_argument("--mat", help="Path to result.mat file from PleioFDR ouput", required=True)
requiredNamed.add_argument("--ref", help="Path to .ref file", required=True)
parser.add_argument("--out", help="Path to file after conversion is done (default: result.mat.csv)")
parser.add_argument("--head", default=5, type=int, help="Number of lines to show (default: 5)")
parser.add_argument("--compress", default=False, action='store_true', help="Compress to .gz archive (default: False)")
args = parser.parse_args()

matfile = args.mat
reffile = args.ref
if args.out is not None:
    outname = args.out
else:
    outname = matfile + '.csv'

# -- Convert result.mat -------------------------

if __name__ == '__main__':
    print('Load {}'.format(matfile))
    sumstats = sio.loadmat(matfile)
    
    print('Load {}'.format(reffile))
    ref = pd.read_csv(reffile, delim_whitespace=True)
    ref['FDR'] = sumstats['fdrmat']

    print('Write {}'.format(outname))
    ref[['CHR', 'SNP', 'BP', 'A1', 'A2', 'FDR']].to_csv(outname, index=False, sep='\t')

    if args.head:
        print(ref.head(args.head))

    if args.compress:
        print('Compress {}'.format(outname))
        os.system('gzip {}'.format(outname))

    print('Done!')

pass
