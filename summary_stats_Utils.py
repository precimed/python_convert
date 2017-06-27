from __future__ import division
import sys, os, re, logging, datetime
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from pkg_resources import parse_version

pdlow = parse_version(pd.__version__) < parse_version('0.17.0')

def read_sumdata(sumFile, snpCol, pCol, kargs):
    '''
    Reading GWAS summary statistics from given file.

    Input:
    ------
    sumFile,    Summary statistics text file, compressed or not
    snpCol,     Field name for SNP identifiers
    pCol,       Field name for P values
    kargs,      Key-words arguments for other information:
                    effCol,     Field name for effect size 
                    ORCol,      Field name for Odds ratio
                    effACol,    Field name for effective allele 
                    othACol,    Field name for other allele 
                    posCol,     Field name for genomic position
                    infoCol,    Field name for imputation quality score
                    NCol,       Field name for sample size
                    sep,        Field separator
    Return:
    -------
    DataFrame,  Field names (if exists) will be standardize according to bellow
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

    '''
    if type(kargs) != dict:
        try:
            kargs = vars(kargs)
        except:
            raise
    if not os.access(sumFile, os.R_OK):
        raise ValueError, 'Unable to read summary file: %s' % (sumFile)
    try:
        if 'sep' in kargs:
            sumDat = pd.read_table(sumFile,sep=kargs['sep'],
                    na_values=[' ', '#N/A','\N','N/A','NA','NULL','NaN', 'nan'])
        else:
            sumDat = pd.read_table(sumFile)
        if sumDat.shape[1] <3: 
            sumDat = pd.read_table(sumFile, sep=' *',
                    na_values=[' ','#N/A','\N','N/A','NA','NULL','NaN', 'nan'])
            if sumDat.shape[1] < 3:
                sumDat = pd.read_table(sumFile, sep='[ +|\t]', engine='python',
                    na_values=[' ', '#N/A','\N','N/A','NA','NULL','NaN', 'nan'])
                if sumDat.shape[1] < 3:
                    raise (ValueError, 
                      "Can't figure out delimiter in %s: tab or space" % (
                          sumFile,))
    except:
        raise
    try:
        sumDat.loc[:,pCol] = sumDat.loc[:,pCol].astype('float') 
        sumDat.rename(columns={snpCol:'SNP', pCol:'P'},inplace=True)
        print sumDat.columns[0]
        misIdx = sumDat.SNP.isnull() | sumDat.P.isnull()
        for k, v in kargs.iteritems():
            if v== None:
                continue
            if k == 'effCol': 
                sumDat.loc[:,v] = sumDat.loc[:, v].astype('float') 
                if v != 'Beta' and 'Beta' in sumDat.columns:
                    sumDat.drop('Beta', axis=1, inplace=True)
                sumDat.rename(columns={v:'Beta'},inplace=True)
                misIdx = misIdx | sumDat.Beta.isnull()
            if k == 'ORCol': 
                if v != 'OR' and 'OR' in sumDat.columns:
                    sumDat.drop('OR', axis=1, inplace=True)
                sumDat.loc[:,v] = sumDat.loc[:, v].astype('float') 
                sumDat.rename(columns={v:'OR'},inplace=True)
                misIdx = misIdx | sumDat.OR.isnull()
            elif k=='infoCol':
                if v != 'INFO' and 'INFO' in sumDat.columns:
                    sumDat.drop('INFO', axis=1, inplace=True)
                sumDat.loc[:,v] = sumDat.loc[:, v].astype('float') 
                sumDat.rename(columns={v:'INFO'},inplace=True)
            elif k=='effACol': 
                if v != 'A1' and 'A1' in sumDat.columns:
                    sumDat.drop('A1', axis=1, inplace=True)
                sumDat.loc[:,v] = sumDat.loc[:, v].str.upper()
                sumDat.rename(columns={v:'A1'},inplace=True)
            elif k=='othACol': 
                if v != 'A2' and 'A2' in sumDat.columns:
                    sumDat.drop('A2', axis=1, inplace=True)
                sumDat.loc[:,v] = sumDat.loc[:, v].str.upper()
                sumDat.rename(columns={v:'A2'},inplace=True)
            elif k=='posCol': 
                if v != 'POS' and 'POS' in sumDat.columns:
                    sumDat.drop('POS', axis=1, inplace=True)
                sumDat.loc[:,v].fillna(-9, inplace=True)
                sumDat.loc[:,v] = sumDat.loc[:, v].astype('int')
                sumDat.rename(columns={v:'POS'},inplace=True)
            elif k=='NCol': 
                if v != 'N' and 'N' in sumDat.columns:
                    sumDat.drop('N', axis=1, inplace=True)
                sumDat.loc[:,v].fillna(-9, inplace=True)
                sumDat.loc[:,v] = sumDat.loc[:, v].astype('float').astype('int')
                sumDat.rename(columns={v:'N'},inplace=True)
            elif k=='chrCol':
                if v != 'CHR' and 'CHR' in sumDat.columns:
                    sumDat.drop('CHR', axis=1, inplace=True)
                sumDat.loc[:,v].fillna(-9, inplace=True)
                sumDat.loc[:,v] = format_chr(sumDat.loc[:,v])
                sumDat.rename(columns={v:'CHR'},inplace=True)
            else:
                pass # leave other names as it is
    except KeyError:
        raise
    except:
        raise #ValueError, 'effect and/or p value is not numeric'
    return sumDat.loc[misIdx==False,:]

def format_chr(chrvec):
    '''
    Reformat chromosome names.

    Input:
    ------
    Vector of chromosome IDs

    Output:
    -------
    Vector of cleaned chromosome IDs

    Note:
    * Remove "chr/Chr/CHR/MT/mt" strings in the name
    * Change chrX to 23, ChrY to 24, PAR to 25, MT to 26
    * (as in plink, https://www.cog-genomics.org/plink/1.9/input#allow_extra_chr)
    '''
    try:
        tmpchrvec = chrvec.astype('str')
        tmpchrvec[tmpchrvec=='X'] = '23'
        tmpchrvec[tmpchrvec=='x'] = '23'
        tmpchrvec[tmpchrvec=='Y'] = '24'
        tmpchrvec[tmpchrvec=='y'] = '24'
        tmpchrvec[tmpchrvec=='PAR'] = '25'
        tmpchrvec[tmpchrvec=='par'] = '25'
        tmpchrvec[tmpchrvec=='M'] = '26'
        tmpchrvec[tmpchrvec=='m'] = '26'
        tmpchrvec[tmpchrvec=='MT'] = '26'
        tmpchrvec[tmpchrvec=='mt'] = '26'
        tmpchrvec = tmpchrvec.str.replace('[chrCHR]', '', case=False)
        # TO-DO: Bellow is anoying
        tmpchrvec[tmpchrvec=='NA'] = '-9'
        tmpchrvec[tmpchrvec.isnull()] = '-9'
        tmpchrvec[tmpchrvec=='nan'] = '-9'
        tmpchrvec[tmpchrvec==' '] = '-9'
        tmpchrvec = tmpchrvec.astype('float').astype('int')
        return tmpchrvec
    except:
        raise

def deduplcate_sum(sumDat, pCol, keys):
    '''
    Remove duplicated rows in give data.

    Input:
    ------
    sumDat,     DataFrame of summary statistics
    pCol,       Field names for p value
    keys,       List of column names used to find duplicated rows

    Return:
    ------
    cleanDat,   DataFrame without duplicated rows.
    dupDat,     DataFrame with only duplicated rows.

    Note:
    -----
    * Only the row with minimal p value will be added to cleanDat.
    * output data will be sorted by give columns

    '''
    if type(keys) == str:
        keys = [keys]
    if pCol not in sumDat.columns:
        raise ValueError, '%s not in columns names' % (pCol,)
    for k in keys:
        if k not in sumDat.columns:
            raise ValueError, '%s not in columns names' % (k,)
    dup_idx = sumDat.duplicated(subset=keys, keep=False)
    cleanDat = sumDat.loc[dup_idx==False,:]
    dup = sumDat.loc[dup_idx,:]
    if dup.shape[0] > 2:
        for k in keys:
            if (np.sum(dup.loc[:,k].isnull())== dup.shape[0]) or \
                    (np.sum(dup.loc[:,k]==-9)== dup.shape[0]) :
                return cleanDat.sort_values(by=list(keys)), dup
        minPDat = sumDat.loc[dup.groupby(by=list(keys))[pCol].idxmin(),:]
        cleanDat = cleanDat.append(minPDat)
    return cleanDat.sort_values(by=list(keys)), dup.sort_values(by=list(keys))

def map_snps(dat1, dat2, keys, suffix, clean=True): 
    '''
    Find data item in/not first DataFrame(dat1) also in second DataFrame(dat2).

    Input:
    ------
    dat1,       First DataFrame.
    dat2,       Second DataFrame.
    keys,       List of column names used to find overlapping rows

    Return:
    ------
    mDat,       DataFrame with common items
    missDat,    DataFrame with only item only existing in dat1.

    Note:
    -----
    * Columns must exists in both DataFrame
    * Assume both DataFrame have already been deduplicated by keys
    '''
    if type(keys) == str:
        keys = [keys]
    for k in keys:
        if k not in dat1.columns:
            raise ValueError, '%s not in columns names of dat1' % (k,)
        if k not in dat2.columns:
            raise ValueError, '%s not in columns names of dat2' % (k,)
    d2names = [k for k in dat2.columns if k not in dat1.columns]
    if d2names == []:
        d2names = [k+'_'+suffix for k in dat2.columns if k not in keys]
    mDat = pd.merge(dat1, dat2, on=keys, how='left', sort=False,
            suffixes=('','_'+suffix))
    missIdx = mDat.loc[:, [k for k in d2names]].isnull().all(axis=1)
    if clean:
        return mDat.loc[missIdx==False, :], mDat.loc[missIdx==True,:]
    else:
        return mDat, mDat.loc[missIdx==True,:]

def flip_snps(sumDat, suffix): 
    '''
    Flip strand for merged sumData.

    Input:
    ------
    sumDat,     Summary statistics DataFrame.
    suffix,     Suffix for the common column names.

    Return:
    ------
    mDat,       DataFrame with common items
    missDat,    DataFrame with only item only existing in dat1.

    Note:
    -----
    * A1, A2 and A1_suffix and A2_suffix must exists in sumDat
    * original A1 and A2 will be aligned with A1_suffix and A2_suffix
    ** TO-DO
    **   may be risky but the pool-man's hope
    '''
    tmpNames = sumDat.columns
    refA1Col = 'A1_'+suffix
    refA2Col = 'A2_'+suffix
    if ('A1' not in tmpNames) or ('A2' not in tmpNames) or \
            (refA1Col not in tmpNames) or (refA2Col not in tmpNames):
        raise (RuntimeError, 'Cant check strand without knowing A1, A2')
    nonflipIdx1 = (sumDat.A1 == sumDat.loc[:,refA1Col]) & (sumDat.A2 ==
            sumDat.loc[:,refA2Col])
    nonflipIdx2 = (sumDat.A1 == sumDat.loc[:,refA2Col]) & (sumDat.A2 ==
            sumDat.loc[:,refA1Col])
    flipIdx = (nonflipIdx1==False) & (nonflipIdx2==False)
    sumDat.loc[flipIdx, 'A1'] = sumDat.A1[flipIdx].map({'A':'T', 'T':'A',
        'C':'G', 'G':'C'})
    sumDat.loc[flipIdx, 'A2'] = sumDat.A2[flipIdx].map({'A':'T', 'T':'A',
        'C':'G', 'G':'C'})
    flipedIdx1 = (sumDat.A1 == sumDat.loc[:,refA1Col]) & (sumDat.A2 ==
            sumDat.loc[:,refA2Col])
    flipedIdx2 = (sumDat.A1 == sumDat.loc[:,refA2Col]) & (sumDat.A2 ==
            sumDat.loc[:,refA1Col])
    flipedIdx = (flipedIdx2==False) & (flipedIdx1==False)
    #assert np.sum(flipedIdx) ==0
    return sumDat

def basic_QC_P(sumDat, outdir, pCol='P', logger=None):
    '''
    Perform qc based on P value information in the summary data.

    Input:
    -----
    sumDat,     Summary statistics DataFrame
    outdir,     Directory of saving intermediate result
    pCol,       Field name of p value
    logger,     object logging to log progress (None)
    
    Return:
    ------
    Note:
    -----
    Remove SNPs with p value >1 or <0, or infinite 

    TO-DO:
    ------
    Add more statistical filters

    '''
    if pCol not in sumDat.columns:
        raise (ValueError, '{} not in the data frame.'.format(pCol))
    Idx = ((sumDat.loc[:, pCol] <=1.0) & (sumDat.loc[:, pCol]>=0)).values
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    outfile = os.path.join(outdir, 'Invalid_P_SNPs.txt.gz')
    logger.info ('Filter out SNPs with invalid pvalues:\n-----------------')
    logger.info ('{} SNPs with p value invalid'.format(np.sum(Idx==False)))
    logger.info ('\n save invalid SNPs to "{}"'.format(outfile))
    logger.info('\n-------------------------------------------')
    tmpD = sumDat.loc[Idx==False,:]
    tmpD.to_csv(outfile, na_rep='NA', compression='gzip', index=False, sep='\t')
    return (sumDat.loc[Idx==True])

def basic_QC_Info(sumDat, outdir, infoCol='INFO', 
        infoTh=0.5, pCol='P',logger=None):
    '''
    Perform qc based on P value information in the summary data.

    Input:
    -----
    sumDat,     Summary statistics DataFrame
    outdir,     Directory of saving intermediate result
    infoCol,    Field name of Imputation info
    infoTh,     Threshold of good impuation
    pCol,       Field name of P value
    logger,     object logging to log progress (None)
    
    Return:
    ------
    Note:
    -----
    Remove SNPs with INFO value < InfoTh or infinite 

    TO-DO:
    ------
    Add more statistical filters

    '''
    if infoCol not in sumDat.columns:
        raise (ValueError, '{} not in the data frame.'.format(infoCol))
    #Idx = ((sumDat.loc[:,infoCol] <= 1.05) & (sumDat.loc[:,infoCol] >= infoTh))
    Idx = (sumDat.loc[:,infoCol] >= infoTh)
    Idx = Idx.values
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    tag = np.empty((sumDat.shape[0],), dtype='|S10'); tag.fill('Clean')
    tag[Idx==False] = 'Removed'
    tmpdf = pd.DataFrame({'Info':sumDat.loc[:,infoCol], 
        '-log10P':-np.log10(sumDat.loc[:,pCol])})
    tmpdf.loc[:, 'Tag'] = tag
    g = sns.FacetGrid(tmpdf, col='Tag', sharex=False, sharey=False)
    g.map(plt.scatter, "Info", "-log10P")
    g.savefig(os.path.join(outdir, 'Low_Info_SNPs.png'))
    g2 = sns.FacetGrid(tmpdf, col='Tag', sharex=False, sharey=False)
    g2.map(plt.hist, "Info")
    g2.savefig(os.path.join(outdir, 'Low_Info_SNPs_hist.png'))
    outfile = os.path.join(outdir, 'Low_INFO_SNPs.txt.gz')
    logger.info ('Filter out SNPs with Imputation info:\n-----------------')
    logger.info ('{} SNPs with low INFO'.format(np.sum(Idx==False)))
    logger.info ('\n save low INFO SNPs to "{}"'.format(outfile))
    logger.info('\n-------------------------------------------')
    tmpD = sumDat.loc[Idx==False,:]
    tmpD.to_csv(outfile, na_rep='NA', compression='gzip', index=False, sep='\t')
    return (sumDat.loc[Idx==True,:])

def basic_QC_Freq(sumDat, outdir, freqCol='Freq', freqTh=0.05, pCol='P',
        logger=None):
    '''
    Perform qc based on frequency of effective allele.

    Input:
    -----
    sumDat,     Summary statistics DataFrame
    outdir,     Directory of saving intermediate result
    freqCol,    Field name of allele frequence
    freqTh,     Frequence threshold
    logger,     object logging to log progress (None)
    
    Return:
    ------
    Note:
    -----
    Remove SNPs with effective allele frequence < freqTh 
    Give up QC if 80% of data failed QC

    TO-DO:
    ------
    Add more statistical filters

    '''
    if freqCol not in sumDat.columns:
        raise (ValueError, '{} not in the data frame.'.format(freqCol))
    freq = sumDat.loc[:,freqCol].values
    Idx = ((freq >= freqTh) & (freq <= 1-freqTh))
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    tmpdf = pd.DataFrame({'H':2*(freq)*(1-freq), 
        '-log10P':-np.log10(sumDat.loc[:,pCol])})
    tag = np.empty((sumDat.shape[0],), dtype='|S10'); tag.fill('Clean')
    tag[Idx==False] = 'Removed'
    tmpdf.loc[:, 'Tag'] = tag
    g = sns.FacetGrid(tmpdf, col='Tag', sharex=False, sharey=False)
    g.map(plt.scatter, "H", "-log10P")
    g.savefig(os.path.join(outdir, 'Maf_filter_SNPs.png'))
    g2 = sns.FacetGrid(tmpdf, col='Tag', sharex=False, sharey=False)
    g2.map(plt.hist, "H")
    g2.savefig(os.path.join(outdir, 'Maf_filter_SNPs_hist.png'))
    outfile = os.path.join(outdir, 'MAF_filtered_SNPs.txt.gz')
    logger.info ('Filter out SNPs with effective allele frequence:\n---------')
    logger.info ('{} SNPs with frequence < {} '.format(np.sum(Idx==False),
        freqTh))
    logger.info ('\n save removed SNPs  to "{}"'.format(outfile))
    logger.info('\n-------------------------------------------')
    if np.sum(Idx==True) <= (sumDat.shape[0]*0.2):
        logger.info('More than 80% of data does not have Freq, deleted column')
        sumDat.loc[:, freqCol] = None
        return (sumDat)
    else:
        tmpD = sumDat.loc[Idx==False,:]
        tmpD.to_csv(outfile, na_rep='NA',compression='gzip',
                index=False, sep='\t')
        return (sumDat.loc[Idx==True,:])

def basic_QC_Freq2(sumDat, outdir, freqACol, freqUCol, freqTh=0.8, 
        NcasCol='', NconCol='', pCol='P', logger=None):
    '''
    Perform qc based on effective allele frequence in Cases and controls.

    Input:
    -----
    sumDat,     Summary statistics DataFrame
    outdir,     Directory of saving intermediate result
    freqACol,   Field name of frequence of effective allele in cases
    freqUCol,   Field name of frequence of effective allele in controls
    freqTh,     Threshold of frequency of effective allele
    NcasCol,    The column name of the number of cases
    NcasCol,    The column name of the number of controls
    logger,     object logging to log progress (None)
    
    Return:
    ------
    Note:
    -----
    Remove SNPs with frequence of effective allele < freqTh in cases and
    conrtols 
    Give up QC when 80% of the data failed QC

    TO-DO:
    ------
    Add more statistical filters

    '''
    if freqACol not in sumDat.columns:
        raise (ValueError, '{} not in the data frame.'.format(freqACol))
    if freqUCol not in sumDat.columns:
        raise (ValueError, '{} not in the data frame.'.format(freqUCol))
    if pCol not in sumDat.columns:
        raise (ValueError, '{} not in the data frame.'.format(pCol))
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    if NcasCol and NconCol:
        freq = np.true_divide(sumDat.loc[:,NcasCol] * sumDat.loc[:,
            freqACol] + sumDat.loc[:, NconCol] * sumDat.loc[:,freqUCol], 
                    sumDat.loc[:, NcasCol]  + sumDat.loc[:, NconCol]) 
    else:
        freq = sumDat.loc[:, freqUCol].values
    Idx = ((freq >= freqTh) & (freq <= (1-freqTh))).values
    tmpdf = pd.DataFrame({'H':2*(freq)*(1-freq), 
        '-log10P':-np.log10(sumDat.loc[:,pCol])})
    tag = np.empty((sumDat.shape[0],), dtype='|S10'); tag.fill('Clean')
    tag[Idx==False] = 'Removed'
    tmpdf.loc[:, 'Tag'] = tag
    logger.info('Filter out SNPs with allele frequences A&U:\n------------')
    logger.info ('{} SNPs with frequences in A < {}'.format(
        np.sum(Idx==False), freqTh))
    outfile = os.path.join(outdir, 'MAF_filtered_SNPs_AU.txt.gz')
    logger.info ('\n save removed  SNPs to "{}"'.format(outfile))
    logger.info('\n-------------------------------------------')
    if np.sum(Idx==True) < (sumDat.shape[0]* 0.2):
        logger.info('More than 80% of data failed QC freq, deleted Column')
        g = sns.FacetGrid(data=tmpdf, col='Tag', sharex=False, sharey=False)
        g.map(plt.scatter, "H", "-log10P")
        g.savefig(os.path.join(outdir, 'Maf_filter_SNPs_AU.png'), 
                width=10, height=10)
        return (sumDat)
    else:
        g1 = sns.FacetGrid(tmpdf, col='Tag', sharex=False, sharey=False)
        g1.map(plt.scatter, "H", "-log10P")
        g1.savefig(os.path.join(outdir, 'Maf_filter_SNPs_AU.png'))
        g2 = sns.FacetGrid(tmpdf, col='Tag', sharex=False, sharey=False)
        g2.map(plt.hist, "H")
        g2.savefig(os.path.join(outdir, 'Maf_filter_SNPs_AU_hist.png'), 
                width=10, height=10)
        tmpD = sumDat.loc[Idx==False,:]
        tmpD.to_csv(outfile, na_rep='NA',compression='gzip',
                index=False,sep='\t')
        return(sumDat.loc[Idx==True,:])

def basic_QC_SNP_only(sumDat, outdir, snpCol='SNP', effACol='A1', 
        othACol='A2', logger=None):
    '''
    Perform qc on selecting SNP only  .

    Input:
    -----
    sumDat,     Summary statistics DataFrame
    outdir,     Directory of saving intermediate result
    snpCol,     SNP ID field name 
    effACol,    effective Allele column
    othACol,    the other Allele column
    logger,     object logging to log progress (None)
    
    Return:
    ------
    Note:
    -----
    Remove Incertion & deletions, CNV
    Assuming SNPs are all with rs number
    '''

    rsIdx = sumDat.loc[:, snpCol].str.match('rs', case=False).values
    A1IIdx = sumDat.loc[:, effACol].str.match('I', case=False).values
    A1DIdx = sumDat.loc[:, effACol].str.match('D', case=False).values
    validIdx = (rsIdx==True) & (A1IIdx==False) & (A1DIdx==False)
    if othACol in sumDat.columns:
        A2IIdx = sumDat.loc[:, othACol].str.match('I', case=False)
        A2DIdx = sumDat.loc[:, othACol].str.match('D', case=False)
        validIdx = (validIdx==True) & (A2IIdx==False) & (A2DIdx == False)
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    outfile = os.path.join(outdir, 'None_SNP.txt.gz')
    logger.info ('Filter out none SNPs with :\n-----------------')
    logger.info ('{} none SNPs variants'.format(np.sum(validIdx==False)))
    logger.info ('\n save none SNPs to "{}"'.format(outfile))
    logger.info('\n-------------------------------------------')
    tmpD = sumDat.loc[validIdx==False,:]
    tmpD.to_csv(outfile, na_rep='NA', compression='gzip', index=False, sep='\t')
    return (sumDat.loc[validIdx==True,:])

def basic_QC_direct(sumDat, outdir, dirCol='DIR', dirMth=1, dirDth=.5, 
        dirPos='+', dirNeg='-', dirMis='?', logger=None):
    '''
    Perform qc based on direction among substudies.

    Input:
    -----
    sumDat,     Summary statistics DataFrame
    outdir,     Directory of saving intermediate result
    dirCol,     Field name of direction
    dirMth,      Threshold of prop of substudies missing direction info
    dirDth,      Threshold of prop of substudies wth discordant direction info
    dirPos,     Symbol for positive direction
    dirNeg,     Symbol for negative direction
    dirMis,     Symbol for missing direction information
    logger,     object logging to log progress (None)
    
    Return:
    ------
    Note:
    -----
    SNPs with too many missing direction removed, i.e. > dirMth 
    SNPs with too much inconsistent direction info removed, i.e., >dirDth

    TO-DO:
    ------
    Add more statistical filters (hetergeniety test?)

    '''
    if dirCol not in sumDat.columns:
        raise (ValueError, '{} not in the data frame.'.format(dirCol))
    nstudy = sumDat.loc[:,dirCol].str.len()
    if nstudy.values[0] <= 6:
        logger.info ('< 6 substudies, skip filter on direction\n')
        return sumDat
    plus = sumDat.loc[:,dirCol].str.count('\\'+dirPos)
    minus = sumDat.loc[:,dirCol].str.count('\\'+dirNeg)
    miss = sumDat.loc[:,dirCol].str.count('\\'+dirMis)
    disRate = np.true_divide(np.min([plus, minus], axis=0), nstudy)
    misRate = np.true_divide(miss, nstudy)
    # be consistent with Ricopili, only consider missing rates
    #Idx = ((misRate <= dirMth) & (disRate <= dirDth))
    Idx = misRate <= dirMth
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    outfile = os.path.join(outdir, 'Invalid_direction_SNPs.txt.gz')
    logger.info ('Filter out SNPs with inconsistent direction:\n------------')
    logger.info ('{} SNPs with inconsistent direction'.format(
        np.sum(Idx==False)))
    logger.info ('\n save invalid SNPs to "{}"'.format(outfile))
    logger.info('\n-------------------------------------------')
    tmpD = sumDat.loc[Idx==False,:]
    tmpD.to_csv(outfile, na_rep='NA', compression='gzip', index=False, sep='\t')
    return (sumDat.loc[Idx==True,:])

def QC_P_Eff(sumDat, outdir, pCol, effCol=None, seCol=None, NCol=None,
        NcasCol=None, NconCol=None, ORCol=None,
        thresh=np.finfo(float), logger=None):
    '''
    Perform qc on the consistency between P and Zscore.

    Input:
    -----
    sumDat,     Summary statistics DataFrame
    outdir,     Directory of saving intermediate result
    effCol,     Column name of effect
    seCol,      Column name of SE
    NCol,       Column name of Sample size
    NcasCol,    Column name of number of cases
    NconCol,    Column name of number of control
    ORCol,      Column name of odds ratio
    thresh,     Threshold of distinguish difference of two floating number
    logger,     object logging to log progress (None)
    
    Return:
    ------
    Note:
    -----
    Remove SNPs having inconsistent effect vs. p
    Check difference between -log10 P and -log10 P(from effect size).
    If se not availbel, sqrt(N) or sqrt(effN) used
    difference of abs(-log10P - (-log10(Peff))) > 1.0, considered as
        inconsistent

    Author: Yunpeng Wang, yunpeng.wng@gmail.com

    '''
    if seCol:
        sevec = sumDat.loc[:, seCol].values
    else:
        if NCol:
            sevec = np.sqrt(sumDat.loc[:, NCol].values)
        elif NcasCol and NconCol:
            effN = 2.0 / ((1.0/sumDat.loc[:, NcasCol].values) +
                    (1.0/sumDat.loc[:,NconCol].values))
            sevec = np.sqrt(effN)
        else:
            sevec = np.ones((sumDat.shape[0],))
    if effCol:
        effvec = sumDat.loc[:, effCol].values
    elif ORCol:
        effvec = np.log(sumDat.loc[:, ORCol].values)
    tmpP = 2*stats.norm.cdf(-np.abs(effvec/sevec))
    valIdx = np.abs(-np.log10(tmpP) - 
            (-np.log10(sumDat.loc[:, pCol].values))) <= thresh
    tag = np.empty((sumDat.shape[0],), dtype='|S10'); tag.fill('Clean')
    tag[valIdx==False] = 'Removed'
    tmpdf = pd.DataFrame({'P':sumDat.loc[:, pCol],'Peff':tmpP, 
        'Tag':tag})
    g = sns.FacetGrid(tmpdf, col='Tag', sharex=False, sharey=False)
    g.map(plt.scatter, "P", "Peff")
    g.savefig(os.path.join(outdir, 'Inconsisten_PorZ_SNP.png'))
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    if np.sum(valIdx==True) < (0.2 * sumDat.shape[0]):
        logger.error('More that 80% data fails QC Eff')
        raise ('Most SNPs failed PvsEff QC')
    else:
        outfile = os.path.join(outdir, 'Inconsisten_PorZ_SNP.txt.gz')
        logger.info ('Filter out SNPs with inconsistent p/z :\n-----------')
        logger.info ('{} SNPs invalid p/z'.format(np.sum(valIdx==False)))
        logger.info ('\n save SNPs with inconsistent P/Z to "{}"'.format(
            outfile))
        logger.info('\n-------------------------------------------')
        tmpD = sumDat.loc[valIdx==False,:]
        tmpD.to_csv(outfile, na_rep='NA', compression='gzip', 
                index=False, sep='\t')
        return (sumDat.loc[valIdx==True,:])
    
def basic_QC_ambiA(sumDat, outdir, effACol='A1', othACol='A2', logger=None):
    '''
    Perform qc on selecting SNP only  .

    Input:
    -----
    sumDat,     Summary statistics DataFrame
    outdir,     Directory of saving intermediate result
    effACol,    effective Allele column
    othACol,    the other Allele column
    logger,     object logging to log progress (None)
    
    Return:
    ------
    Note:
    -----
    Remove SNPs having ambiguouse allele coding AT, CG
    Assuming SNPs are all with rs number
    '''

    A1AIdx = sumDat.loc[:, effACol].str.match('A', case=False)
    A1GIdx = sumDat.loc[:, effACol].str.match('G', case=False)
    A1CIdx = sumDat.loc[:, effACol].str.match('C', case=False)
    A1TIdx = sumDat.loc[:, effACol].str.match('T', case=False)

    A2AIdx = sumDat.loc[:, othACol].str.match('A', case=False)
    A2GIdx = sumDat.loc[:, othACol].str.match('G', case=False)
    A2CIdx = sumDat.loc[:, othACol].str.match('C', case=False)
    A2TIdx = sumDat.loc[:, othACol].str.match('T', case=False)

    ambiIdx = (A1AIdx==True) & (A2TIdx==True)
    ambiIdx = (ambiIdx==True ) | ((A1TIdx==True) & (A2AIdx==True))
    ambiIdx = (ambiIdx==True ) | ((A1CIdx==True) & (A2GIdx==True))
    ambiIdx = (ambiIdx==True ) | ((A1GIdx==True) & (A2CIdx==True))
    ambiIdx = ambiIdx.values
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    outfile = os.path.join(outdir, 'Ambiguous_allele_SNP.txt.gz')
    logger.info ('Filter out SNPs with ambiguous allele coding:\n-----------')
    logger.info ('{} SNPs with AT/TA/CG/GC'.format(np.sum(ambiIdx==True)))
    logger.info ('\n save ambiguous SNPs to "{}"'.format(outfile))
    logger.info('\n-------------------------------------------')
    tmpD = sumDat.loc[ambiIdx==True,:]
    tmpD.to_csv(outfile, na_rep='NA', compression='gzip', index=False, sep='\t')
    return (sumDat.loc[ambiIdx==False,:])

def deduplicate_bycol(dat,  keys):
    '''
    Remove duplicated rows in give data.

    Input:
    ------
    dat,        DataFrame 
    keys,       List of column names used to find duplicated rows

    Return:
    ------
    cleanDat,   DataFrame without duplicated rows.
    dupDat,     DataFrame with only duplicated rows.

    Note:
        Keep the first.


    '''
    if type(keys) == str:
        keys = [keys]
    for k in keys:
        if k not in dat.columns:
            raise ValueError, '%s not in columns names' % (k,)
    dupIdx = dat.duplicated(subset=keys)
    cleanDat = dat.loc[dupIdx==False,:]
    dupDat = dat.loc[dupIdx,:]
    return cleanDat, dupDat

