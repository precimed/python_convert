import pandas as pd 
import numpy as np
import os, sys, logging

def read_bim(bimFile, logger=None, sep='\t'):
    '''
    Read PLINK bim file.

    Input:
    ------
    bimFile,    PLINK bim file path
    logger,     python logger for process information
    sep,        separator of bim file

    Return:
    ------
    bimDat,     DataFrame with SNP information

    Note:
    -----
    * Adding column names for convienience:
        CHR SNP GP  POS A1  A2
    * Change ChrX->23, ChrY->24 and ChrM->25
    '''
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    if not os.access(bimFile, os.R_OK):
        logger.error('Unable to read {}'.format(bimFile))
        raise (ValueError, 'Unable to read {}'.format(bimFile))
    bimDat = pd.read_csv(bimFile, sep=sep, header=None, 
            names=['CHR', 'SNP', 'GP', 'POS', 'A1', 'A2'])
    bimDat.loc[:,'CHR'] = bimDat.loc[:, 'CHR'].astype('|S5')
    bimDat.loc[bimDat.loc[:, 'CHR']=='X', 'CHR'] = '23'
    bimDat.loc[bimDat.loc[:, 'CHR']=='Y', 'CHR'] = '24'
    bimDat.loc[bimDat.loc[:, 'CHR']=='M', 'CHR'] = '25'
    bimDat.loc[:, 'CHR'] = bimDat.loc[:, 'CHR'].astype('float').astype('int')
    bimDat.loc[:, 'POS'] = bimDat.loc[:, 'POS'].astype('int')
    logger.info('Read {} SNPs from {}'.format(bimDat.shape[0], bimFile))
    logger.info('Columns: CHR, SNP, GP, POS, A1, A2 were used')
    return (bimDat)

def deduplicate_bim(bimDat, outdir, logger=None):
    '''
    Check if PLINK bim data has duplicate SNP by position.

    Input:
    ------
    bimDat,     DataFrame with PLINK bim data
    outdir,     Output directory for intemediate files
    logger,     python logger for process information

    Return:
    -------
    dupIdx,     Indictor Series for SNPs that should be removed

    Note:
    -----
    * Save all duplicated SNPs into a gziped file.
    * For duplicated SNPs:
    *   SNPs with rs-number kept, i.e., A1 is (A,T,C,G) and A2 is (A,T,C,G) and
    *   SNP ID starts with 'rs'
        Otherwise,
        the first of the duplicates were kept.
    * Also save a text file containing SNP IDs that should be removed by PLINK
    Warning:
    * Extreme slow for large dataset. So better do it once and update the 
    *   corresponding PLINK bed/bim/fam. 
    '''
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    dupIdx = bimDat.duplicated(subset=['CHR', 'POS'], keep=False)
    ndup = np.sum(dupIdx)
    if ndup > 0:
        outfile = os.path.join(outdir, 'Duplicated_SNPs_by_POS.txt.gz')
        logger.warn('Bim file has {} duplicated items by genomic position'.format(ndup))
        logger.warn('Save all duplicated SNPs in Bim to {}'.format(outfile))
        dupDat = bimDat.loc[dupIdx==True,:]
        dupDat.to_csv(outfile, index=False, compression='gzip', na_rep='NA',
                sep='\t')
        grouped = dupDat.groupby(by=['CHR','POS'], sort=False)
        for name, x in grouped:
            rsIdx = x.loc[:,'SNP'].str.startswith('rs')
            A1Idx = x.loc[:,'A1'].str.contains('[A|T|C|G]')
            A2Idx = x.loc[:,'A2'].str.contains('[A|T|C|G]')
            Idx = rsIdx & A1Idx & A2Idx
            if np.sum(Idx) != 1:
                Idx.iloc[0] = False
            dupIdx.values[Idx.index] = Idx.values 
        outfile = os.path.join(outdir, 'Duplicated_SNPs_by_POS_excluded.txt.gz')
        logger.warn('Save removed duplicated SNPs in Bim to {}'.format(outfile))
        tmpDat = dupDat.loc[dupIdx==True]
        tmpDat.to_csv(outfile, index=False, compression='gzip', na_rep='NA',
                sep='\t', columns=['SNP'])
    return (dupIdx) 
