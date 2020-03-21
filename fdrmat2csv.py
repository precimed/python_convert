# Usage:
# python fdrmat2csv.py result.mat /space/syn03/1/data/GWAS/SUMSTAT/misc/9279485_ref.bim

import pandas as pd
import scipy.io as sio
import sys
import numpy as np
if __name__ == '__main__':
    matfile = sys.argv[1]
    reffile = sys.argv[2]
    print('Load {}'.format(matfile))
    sumstats = sio.loadmat(matfile)
    print('Load {}'.format(reffile))
    ref=pd.read_table(reffile, delim_whitespace=True)
    ref['FDR']=sumstats['fdrmat']
    print('Write {}'.format(matfile + '.csv'))
    ref[['CHR', 'SNP', 'BP', 'A1', 'A2', 'FDR']].to_csv(matfile + '.csv', index=False, sep='\t')
    print('Done')
