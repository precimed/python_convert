import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.io as sio
import os, sys, argparse, time, logging, getpass
import matplotlib.pyplot as plt
from summary_stats_Utils import *

def read_sum_dat(sumFile, logger, kargs):
    '''
    Read give summary statistics.

    Input:
    sumFile,    Path of summary file.
    logger,     Python logger for process information.
    kargs,      namespace object for options.

    Return:
    -------
    sumDat,     DataFrame of summary dataset.

    Note:
    -----
    1. Field names (if exists) will be standardize according to bellow
                        effCol,     Beta
                        ORCol,      OR
                        effACol,    A1
                        othACol,    A2
                        posCol,     POS
                        infoCol,    INFO
                        NCol,       N
                And, Chromosome names will be standardized
                        Removing'CHR', 'Chr', etc --> integer
                        recode chrX--> 23
                        recode chrY--> 24
                        recode chrM--> 25
    2. SNPs with invalid p values removed, i.e., >1 or < 0, or NAN
    3. Duplicated SNPs removed
    '''

    if not os.access(sumFile, os.R_OK):
        raise ValueError("Can't read summary stats file: {}".format(sumFile))
    logger.info('*** Loading summary stats ***')
    logger.info('Read summary data from {}'.format(sumFile))
    sumDat = read_sumdata(sumFile, kargs.snpCol, kargs.pCol, kargs)
    logger.info('......')
    logger.info('Read {} SNPs'.format(sumDat.shape[0]))
    colnames = ['SNP', 'P', 'A1', 'CHR', 'POS', 'Beta', 'A2']
    if 'P' not in sumDat.columns:
        raise RuntimeError('No P value provided')
    if 'SNP' not in sumDat.columns:
        raise RuntimeError('No SNP ID provided') 
    if not kargs.effACol:
        warnings.warn('No effective Allele provided') 
        logger.warning('No effective Allele provided') 
        colnames.remove('A1')
    if not kargs.othACol:
        warnings.warn( "No other Allele information provided") 
        logger.warn('No effective Allele provided') 
        colnames.remove('A2')
    if not kargs.effCol:
        if not kargs.orCol:
            colnames.remove('Beta')
            logger.warn('Directionality is not checked')
        else:
            sumDat.loc[:, 'Beta'] = np.log(sumDat.loc[:, 'OR'])
            sumDat.drop('OR', axis=1, inplace=True)
    if (not kargs.effACol) and (not kargs.othACol):
        logger.warn('Directionality is not checked')
        colnames.remove('Beta')
        sumDat.drop('Beta', axis=1, inplace=True)
    if ((not kargs.posCol) or (not kargs.chrCol)) and (not kargs.chrPosCol):
        logger.info('Using SNP ID only for align Summary data to reference')
        colnames.remove('POS')
        colnames.remove('CHR')
        keys = ['SNP']
    elif kargs.forceID:
        keys = ['SNP']
    else:
        keys = ['CHR', 'POS']
    logger.info('Reading Summary stats done\n')
    logger.info('**** check P values *****')
    sumDat = basic_QC_P(sumDat, kargs.O, 'P', logger)
    logger.info('**** END check P values *****')
    logger.info('**** check duplicated SNPs *****')
    sumDat, dup = deduplcate_sum(sumDat, 'P', keys)
    if dup.shape[0] > 0:
        dupFile = os.path.join(kargs.O, 'Duplicated_snps.gz')
        logger.warning('There are {} duplicated SNPs in {}'.format(
            dup.shape[0], sumFile))
        logger.warning('\t The SNP with minimum p value included')
        logger.warning('see all duplicated SNPs in {}'.format(dupFile))
        dup.to_csv(dupFile, index=False, na_rep='NA', compression='gzip',
                sep='\t')
        logger.info('**** END check duplicated SNPs *****')
    sumDat = sumDat.loc[:, colnames] 
    logger.info('\n')

    logger.info('**** save standardized summary stats *****')
    sumDatFile = os.path.join(kargs.O, 'Standardized_snps.gz')
    sumDat.to_csv(sumDatFile, index=False, na_rep='NA', compression='gzip',
                sep='\t')
    logger.info('Save summary data to {}'.format(sumDatFile))
    return sumDat

def read_ref_dat(refFile, logger):
    '''
    Read in-house reference dataset.

    Input:
    ------
    refFile,    Path of reference file:
                CHR, SNP, GP, BP, A1, A2, complementA1, complementA2
    logger,     Python logger for process information

    Return:
    ------
    refDat,     DataFrame of reference dataset.

    Example
    -------
    * https://s3-eu-west-1.amazonaws.com/biostat/2558411_ref.bim.gz
    * https://s3-eu-west-1.amazonaws.com/biostat/9279485_ref.bim.gz
    '''

    if not os.access(refFile, os.R_OK):
        raise ValueError("Can't read reference file: {}".format(refFile))
    logger.info('*** Loading reference data ***')
    refDat = pd.read_table(refFile)
    refDat.rename(columns={'BP':'POS', 'A1':'refA1', 'A2':'refA2'}, 
            inplace=True)
    logger.info('Read reference data from {}'.format(refFile))
    logger.info('Read {} SNPs from reference data'.format(refDat.shape[0]))
    print ('*** Using reference with {} SNPs ***'.format(refDat.shape[0]))
    logger.info('Reading reference data done\n')
    logger.info('\n')
    return(refDat)

def _qq(pvec, ax):
    '''
    Making basic QQ plots of pvalues.

    '''
    pvec = pvec[np.isfinite(pvec)]
    pvec[pvec < 1e-20] = 1e-20
    logpSort =  -np.log10(np.sort(pvec))
    n = logpSort.shape[0]
    logpTheo = -np.log10(np.cumsum(np.repeat(1.0/n, n)))
    ax.scatter(logpTheo, logpSort)
    x = np.linspace(*ax.get_xlim())
    ax.plot(x, x)
    plt.xlabel('Theorectial -log10 (P)')
    plt.ylabel('Observed -log10 (P)')

def summarize_merge(sumDat, mDat, misDat, outdir, logger):
    '''
    Making QQ plot of original dataset, converted and missed.

    Input:
    ------
    sumDat,     DataFrame of Original summary stats
    mDat,       DataFrame of Converted summary data
    misDat,     DataFrame of SNPs in original but not in converted
    outdir,     Where to save figure
    logger,     Python logger for process information

    No return.
    ----------
    TO-DO:
        Making multiple curves in one figure
    '''
    logger.info('\n')
    if sumDat.shape[0] < 10:
        logger.erro('Too few SNPs converted!! N={}'.format(sumDat.shape[0]))
        raise (RuntimeError, 
                'Too few SNPs converted!! N={}'.format(sumDat.shape[0]))
    fig = plt.figure(facecolor='white')
    ax = fig.add_subplot(131)
    _qq(sumDat.loc[:,'P'].values, ax)
    plt.title('Original')
    ax = fig.add_subplot(132)
    _qq(mDat.loc[:,'P'].values, ax)
    plt.title("Converted")
    ax = fig.add_subplot(133)
    _qq(misDat.loc[:,'P'].values, ax)
    plt.title("Missed")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'QQ_convert.png'), format='png')
    plt.close()
    logger.info('Comparing P values in QQ_convert.png')

def check_zscore(zvec, outdir, logger):
    '''
    Check distribution of converted z-score(real not Anders')

    Input:
    ------
    outdir,     Where to save figure
    logger,     Python logger for process information

    No return.
    '''
    logger.info('\n')
    fig = plt.figure(facecolor='white')
    pd.Series(zvec[np.isfinite(zvec)]).hist(bins=100) 
    plt.title('Z-Scores')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'Z_scores.png'), format='png')
    plt.close()
    logger.info('Check converted Z-scores at Z_scores.png')

def align2ref(sumDat, refDat, outdir, logger, kargs):
    '''
    Align given summary Data to in-house reference dataset.

    Input:
    ------
    sumDat,     DataFrame of summary statistics.
    refDat,     DataFrame of in-house reference dataset.
    outdir,     Where to save figure
    logger,     Python logger for process information
    kargs,      NameSpace object of options

    Return:
    -------
    -log10 p values, and z-scores

    Note:
    -----
    1. Ambiguous SNPs removed based on in-house reference dataset.
    2. effect aligned with allele coding of reference
    '''
    if kargs.forceID:
        keys = ['SNP']
    elif ('CHR' not in sumDat.columns) or ('POS' not in sumDat.columns):
        keys = ['SNP']
    else:
        keys = ['CHR', 'POS']
    mDat, misDat1 = map_snps(refDat, sumDat, 'sum', keys, False)
    mDat.to_csv(os.path.join(kargs.O, 'debug_merged.txt.gz'), 
        sep='\t', index=False, na_rep='NA')
    logger.info('*** Align SNPs to reference ***') 
    if misDat1.shape[0] > 0:
        outF = os.path.join(kargs.O, 'SNPs_not_in_sumFile.txt.gz')
        logger.info(
           'There are {} SNPs in reference not in given summary file'.format(
               misDat1.shape[0]))
        logger.info('Details see {}'.format(outF))
        misDat1.to_csv(outF, index=False, sep='\t', compression='gzip',
                na_rep='NA')
    dummy, misDat2 = map_snps(sumDat, refDat, 'ref', keys)
    if misDat2.shape[0] > 0:
        outF = os.path.join(kargs.O, 'SNPs_not_in_refFile.txt.gz')
        logger.info(
            'There are {} SNPs in summary file not in reference'.format(
                misDat2.shape[0]))
        logger.info('Details see {}'.format(outF))
        misDat2.to_csv(outF, index=False, sep='\t', compression='gzip',
                na_rep='NA')
    signvec = np.empty((mDat.shape[0],), dtype='float'); signvec.fill(np.nan)
    ambivec = (((mDat.refA1=='A')&(mDat.refA2=='T')) | 
        ((mDat.refA2=='A')&(mDat.refA1=='T')) |
        ((mDat.refA1=='C')&(mDat.refA2=='G')) |
        ((mDat.refA2=='C')&(mDat.refA1=='G')))
    ambivec = ambivec.values
    logger.info('{} SNPs have ambiguously coded allele in ref'. format(
        np.sum(ambivec)))
    logger.info('Zscores of ambiguously coded SNPs were set to NaN')
    ambDat = mDat.loc[ambivec,:]
    ambDat.to_csv(os.path.join(kargs.O, 'Ambiguous_data.txt.gz'),
            compression='gzip', sep='\t', index=False, na_rep='NA')
    logger.info('Save SNPs with ambiguous allele coding into {}'.format(
        os.path.join(kargs.O, 'Ambiguous_data.txt.gz')))
    logpvec = -np.log10(mDat.loc[:,'P'])
    if 'A1' not in sumDat.columns:
        zvec = signvec.copy()
    else:
        if 'A2' not in sumDat.columns:
            idx1 = ((mDat.A1==mDat.refA1) | (mDat.A1==mDat.A1c)).values
            idx_1 = ((mDat.A1==mDat.refA2) | (mDat.A1==mDat.A2c)).values
        else:
            idx1 = (((mDat.A1==mDat.refA1)&(mDat.A2==mDat.refA2)) | ((mDat.A1==mDat.A1c)&(mDat.A2==mDat.A2c))).values
            idx_1 = (((mDat.A1==mDat.refA2)&(mDat.A2==mDat.refA1)) | ((mDat.A1==mDat.A2c)&(mDat.A2==mDat.A1c))).values
    signvec[idx1] = 1.0; signvec[idx_1] = -1.0; signvec[ambivec] = np.nan
    signvec = signvec * np.sign(mDat.loc[:,'Beta'].values)
    zvec = np.abs(stats.norm.ppf(mDat.loc[:,'P'].values * 0.5)) * signvec
    logger.info('{} SNPs have direction opposite to refference and changed'.format(np.sum(idx_1)))
    mDat.loc[:, 'newZ'] = zvec
    tmpMdat = mDat.loc[idx_1 ,:]
    tmpMdat.to_csv(os.path.join(outdir, 'flip_data.txt.gz'),
        index=False, sep='\t', compression='gzip',na_rep='NA')
    summarize_merge(sumDat, mDat, misDat2, outdir, logger)
    logger.info('\n')
    return(logpvec.values, zvec)

def save2mat(logpvec, zvec, trait, outdir, logger):
    '''
    Save data in Matlab dataset.

    Input:
    -----
    logpvec,    -log10 p value vector
    zvec,       zscore vector
    trait,      Name of phenotype
    outdir,     Where to save dataset
    logger,     Python logger for process information

    No return.
    '''
    outfile = os.path.join(outdir, trait)
    tmpdict = {'logpvec_'+trait.lower():logpvec, 'zvec_'+trait.lower():zvec}
    sio.savemat(outfile, tmpdict, format='5', do_compression=False,
            oned_as='column')
    logger.info('Save converted data to {}'.format(outfile+'.mat'))
    
def convert_sum():
    parser = argparse.ArgumentParser(prog="Preprocess Summary stats",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description='Preprocess summary stats for matlab') 
    parser.add_argument('-F', type=str, help='Summary stats file')
    parser.add_argument('-Ref', type=str, help='Reference file (optional);')
    parser.add_argument('-T', type=str, help='Trait Name')
    parser.add_argument('-O', type=str, help='Output DIR', default=".")
    parser.add_argument('--forceID', action='store_true', default=False,
            help='Force using SNP ID other than position')
    parser.add_argument('--snpCol', type=str, help='SNP ID field', 
            default='SNP')
    parser.add_argument('--pCol', type=str, help='P value field', default='P')
    parser.add_argument('--effACol', type=str, help='Effective allele field', 
            default=None)
    parser.add_argument('--othACol', type=str, help='The other allele field', 
            default=None)
    parser.add_argument('--effCol', type=str, help='Effect size field', 
            default=None)
    parser.add_argument('--orCol', type=str, help='Odds ratio field', 
            default=None)
    parser.add_argument('--posCol', type=str,
            help='Genomic position field',default=None) 
    parser.add_argument('--chrCol', type=str,
            help='Chromosome field',default=None) 
    parser.add_argument('--chrPosCol', type=str,
            help='Field name for a single column with both chromosome and position ' +
                 '(joined by colon, example: "8:103044620")', default=None) 
    args = parser.parse_args()
    if args.O is None: raise ValueError("Output DIR is not provided")
    if args.F is None: raise ValueError("Summary stats file is not provided")
    if not os.access(args.O, os.F_OK):
        os.mkdir(args.O)
    if not os.access(args.F, os.R_OK):
        raise ValueError("Can't read summary stats file: {}".format(args.F))
    if args.Ref and not os.access(args.Ref, os.R_OK):
        raise ValueError("Can't read reference file: {}".format(args.Ref))
    logfile = os.path.join(args.O, 'convert_' + args.T + '.log')
    logger = logging.getLogger()
    logger.addHandler(logging.FileHandler(logfile))
    logger.setLevel(logging.DEBUG)
    sumDat = read_sum_dat(args.F, logger, args)
    if args.Ref:
        refDat = read_ref_dat(args.Ref, logger)
        logpvec, zvec = align2ref(sumDat, refDat, args.O, logger, args)
        check_zscore(zvec, args.O, logger)
        save2mat(logpvec, zvec, args.T, args.O, logger)
    logger.info('\n**********\nFinished at {}'.format(time.ctime()))
    logger.info('Author: {} at {}'.format(getpass.getuser(), time.ctime()))
    

if __name__ == "__main__":
    import time
    import numpy as np
    tsts = time.time()
    convert_sum()
    print
    print ('Finish at {}'.format(time.ctime()))
    ted = time.time()
    print ('Time taken {} mins {} sec'.format((ted-tsts)//60, np.round(ted-tsts) % 60))

