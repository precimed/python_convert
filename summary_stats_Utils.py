import sys, os, re, logging, datetime
import numpy as np
import pandas as pd
from six import iteritems
from pkg_resources import parse_version

pdlow = parse_version(pd.__version__) < parse_version('0.17.0')

def _get_str_list_sign(str_list):
    return np.array([-1 if e[0]=='-' else 1 for e in str_list], dtype=np.int)

def _fix_chromosome_labels(sumDat, v='CHR'):
    sumDat.loc[:,v] = pd.Series([re.sub('[chrCHR]','',
        str(c)) for c in sumDat.loc[:,v]])
    sumDat.loc[sumDat.loc[:, v]=='X', v] = '23'
    sumDat.loc[sumDat.loc[:, v]=='x', v] = '23'
    sumDat.loc[sumDat.loc[:, v]=='Y', v] = '24'
    sumDat.loc[sumDat.loc[:, v]=='y', v] = '24'
    sumDat.loc[sumDat.loc[:, v]=='M', v] = '25'
    sumDat.loc[sumDat.loc[:, v]=='m', v] = '25'
    sumDat.loc[:,v] = sumDat.loc[:, v].astype('float').astype('int')


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
                    chrCol,     Field name for chromosome number
                    posCol,     Field name for genomic position
                    chrPosCol,  Field name for a single column with both chromosome and position
                                (joined by colon, example: "8:103044620")
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
    print(pd.__version__)
    if type(kargs) != dict:
        try:
            kargs = vars(kargs)
        except:
            raise
    if not os.access(sumFile, os.R_OK):
        raise ValueError('Unable to read summary file: %s' % (sumFile))
    try:
        dtype_map = {kargs['effCol']: str} if kargs['effCol'] else {}
        if 'sep' in kargs:
            sumDat = pd.read_table(sumFile, sep=kargs['sep'], dtype=dtype_map)
        else:
            sumDat = pd.read_table(sumFile, dtype=dtype_map)
            print(sumDat.columns)
            print(sumDat.shape)
        if sumDat.shape[1] <2: 
            sumDat = pd.read_table(sumFile, sep=' *', dtype=dtype_map)
            if sumDat.shape[1] < 2:
                sumDat = pd.read_table(sumFile, sep='[ +|\t]', engine='python', dtype=dtype_map)
                if sumDat.shape[1] < 2:
                    raise (ValueError(
                      "Can't figure out delimiter in %s: tab or space" % (
                          sumFile,)))
    except:
        raise
    try:
        sumDat.loc[:,pCol] = sumDat.loc[:,pCol].astype('float') 
        sumDat.rename(columns={snpCol:'SNP', pCol:'P'},inplace=True)
        if kargs['effCol']:
            # effect column has str type
            # -1 if effect starts with '-' else 1
            sumDat.loc[:,'SIGN'] = _get_str_list_sign(sumDat[kargs['effCol']])
        elif kargs['orCol'] :
            # effect column has np.float type
            # 1 if effect >=1 else -1
            if (sumDat[kargs['orCol']] < 0).any():
                raise ValueError("{} column contains negative values; consider using --effCol instead of --orCol".format(kargs['orCol']))
            effect_sign = np.sign(sumDat.loc[:, kargs['orCol']] - 1)
            effect_sign[effect_sign == 0] = 1             # ToBe discussed --- how do we interpred OR == 1.0
            sumDat.loc[:,'SIGN'] = effect_sign
        for k, v in iteritems(kargs):
            print('{} {}'.format(k, v))
            if v== None:
                continue
            if k == 'effCol': 
                sumDat.loc[:,v] = sumDat.loc[:, v].astype('float') 
                sumDat.rename(columns={v:'Beta'},inplace=True)
            if k == 'orCol': 
                sumDat.loc[:,v] = sumDat.loc[:, v].astype('float') 
                sumDat.rename(columns={v:'OR'},inplace=True)
            elif k=='infoCol':
                sumDat.loc[:,v] = sumDat.loc[:, v].astype('float') 
                sumDat.rename(columns={v:'INFO'},inplace=True)
            elif k=='effACol': 
                sumDat.loc[:,v] = sumDat.loc[:, v].str.upper()
                sumDat.rename(columns={v:'A1'},inplace=True)
            elif k=='othACol': 
                sumDat.loc[:,v] = sumDat.loc[:, v].str.upper()
                sumDat.rename(columns={v:'A2'},inplace=True)
            elif k=='posCol': 
                sumDat.loc[:,v] = sumDat.loc[:, v].astype('int')
                sumDat.rename(columns={v:'POS'},inplace=True)
            elif k=='NCol': 
                sumDat.loc[:,v] = sumDat.loc[:, v].astype('int')
                sumDat.rename(columns={v:'N'},inplace=True)
            elif k=='chrCol':
                sumDat.rename(columns={v:'CHR'},inplace=True)
                _fix_chromosome_labels(sumDat)
            elif k=='chrPosCol':
                 sumDat['CHR'], sumDat['POS'] = zip(*sumDat[v].apply(lambda x: x.split(':', 1)))
                 sumDat.loc[:,'POS'] = sumDat.loc[:,'POS'].apply(lambda x: x if x.isdigit() else '-1').astype('int')
                 _fix_chromosome_labels(sumDat)
            else:
                pass # leave other names as it is
    except KeyError:
        raise
    except:
        raise #ValueError, 'effect and/or p value is not numeric'
    return sumDat

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
        raise ValueError('%s not in columns names' % (pCol,))
    for k in keys:
        if k not in sumDat.columns:
            raise ValueError('%s not in columns names' % (k,))
    dup_idx = sumDat.duplicated(subset=keys, keep=False)
    cleanDat = sumDat.loc[dup_idx==False,:]
    dup = sumDat.loc[dup_idx,:]
    if dup.shape[0] > 2:
        minPDat = sumDat.loc[dup.groupby(by=list(keys))[pCol].idxmin(),:]
        cleanDat = cleanDat.append(minPDat)
    return cleanDat.sort_values(by=list(keys)), dup.sort_values(by=list(keys))

def map_snps(dat1, dat2, rsufix, keys, clean=True): 
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
            raise ValueError('%s not in columns names of dat1' % (k,))
        if k not in dat2.columns:
            raise ValueError('%s not in columns names of dat2' % (k,))
    tmpD1 = dat1.set_index(keys = list(keys), drop=False)
    tmpD2 = dat2.set_index(keys = list(keys), drop=False)
    mDat = tmpD1.join(tmpD2, rsuffix=rsufix, how='left', sort=False)
    miss_idx = mDat.loc[:, [k+rsufix for k in keys]].isnull().all(axis=1)
    if clean:
        return mDat.loc[miss_idx==False, :], mDat.loc[miss_idx==True,:]
    else:
        return mDat, mDat.loc[miss_idx==True,:]

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
        raise (ValueError('{} not in the data frame.'.format(pCol)))
    Idx = ((sumDat.loc[:, pCol] <=1.0) & (sumDat.loc[:, pCol]>=0))
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

def basic_QC_Info(sumDat, outdir, infoCol='INFO', infoTh=0.5, logger=None):
    '''
    Perform qc based on P value information in the summary data.

    Input:
    -----
    sumDat,     Summary statistics DataFrame
    outdir,     Directory of saving intermediate result
    infoCol,    Field name of Imputation info
    infoTh,     Threshold of good impuation
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
        raise (ValueError('{} not in the data frame.'.format(infoCol)))
    Idx = ((sumDat.loc[:,infoCol] <=1.05) & (sumDat.loc[:,infoCol]>=infoTh))
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    outfile = os.path.join(outdir, 'Low_INFO_SNPs.txt.gz')
    logger.info ('Filter out SNPs with Imputation info:\n-----------------')
    logger.info ('{} SNPs with low INFO'.format(np.sum(Idx==False)))
    logger.info ('\n save low INFO SNPs to "{}"'.format(outfile))
    logger.info('\n-------------------------------------------')
    tmpD = sumDat.loc[Idx==False,:]
    tmpD.to_csv(outfile, na_rep='NA', compression='gzip', index=False, sep='\t')
    return (sumDat.loc[Idx==True])

def basic_QC_Freq(sumDat, outdir, freqCol='Freq', freqTh=0.05, logger=None):
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

    TO-DO:
    ------
    Add more statistical filters

    '''
    if freqCol not in sumDat.columns:
        raise (ValueError('{} not in the data frame.'.format(freqCol)))
    Idx = ((sumDat.loc[:,freqCol] >= freqTh) & 
            (sumDat.loc[:,freqCol] <= 1-freqTh))
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    outfile = os.path.join(outdir, 'MAF_filtered_SNPs.txt.gz')
    logger.info ('Filter out SNPs with effective allele frequence:\n-----------------')
    logger.info ('{} SNPs with frequence < {} '.format(np.sum(Idx==False),
        freqTh))
    logger.info ('\n save removed SNPs  to "{}"'.format(outfile))
    logger.info('\n-------------------------------------------')
    tmpD = sumDat.loc[Idx==False,:]
    tmpD.to_csv(outfile, na_rep='NA', compression='gzip', index=False, sep='\t')
    return (sumDat.loc[Idx==True])

def basic_QC_Freq2(sumDat, outdir, freqACol, freqUCol, freqTh=0.8, logger=None):
    '''
    Perform qc based on effective allele frequence in Cases and controls.

    Input:
    -----
    sumDat,     Summary statistics DataFrame
    outdir,     Directory of saving intermediate result
    freqACol,   Field name of frequence of effective allele in cases
    freqUCol,   Field name of frequence of effective allele in controls
    freqTh,     Threshold of frequency of effective allele
    logger,     object logging to log progress (None)
    
    Return:
    ------
    Note:
    -----
    Remove SNPs with frequence of effective allele < freqTh in cases and
    conrtols 

    TO-DO:
    ------
    Add more statistical filters

    '''
    if freqACol not in sumDat.columns:
        raise (ValueError('{} not in the data frame.'.format(freqACol)))
    if freqUCol not in sumDat.columns:
        raise (ValueError('{} not in the data frame.'.format(freqUCol)))
    IdxA = ((sumDat.loc[:, freqACol] >= freqTh) & 
             (sumDat.loc[:, freqACol] <= (1-freqTh)))
    IdxU = ((sumDat.loc[:, freqUCol] >= freqTh) & 
            (sumDat.loc[:, freqUCol] <= (1-freqTh)))
    Idx = (IdxA & IdxU)
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    outfile = os.path.join(outdir, 'MAF_filtered_SNPs_AU.txt.gz')
    logger.info ('Filter out SNPs with allele frequences A&U:\n-----------------')
    logger.info ('{} SNPs with frequences in A < {}'.format(np.sum(IdxA==False), freqTh))
    logger.info ('{} SNPs with frequences in U < {}'.format(np.sum(IdxU==False), freqTh))
    logger.info ('\n save removed  SNPs to "{}"'.format(outfile))
    logger.info('\n-------------------------------------------')
    tmpD = sumDat.loc[Idx==False,:]
    tmpD.to_csv(outfile, na_rep='NA', compression='gzip', index=False, sep='\t')
    return(sumDat.loc[Idx==True])

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
        raise (ValueError('{} not in the data frame.'.format(dirCol)))
    plus = sumDat.loc[:,dirCol].str.count('\\'+dirPos)
    minus = sumDat.loc[:,dirCol].str.count('\\'+dirNeg)
    miss = sumDat.loc[:,dirCol].str.count('\\'+dirMis)
    nstudy = sumDat.loc[:,dirCol].str.len()
    disRate = np.true_divide(np.min([plus, minus], axis=0), nstudy)
    misRate = np.true_divide(miss, nstudy)
    Idx = ((misRate <= dirMth) & (disRate <= dirDth))
    if not logger:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler())
    outfile = os.path.join(outdir, 'Invalid_direction_SNPs.txt.gz')
    logger.info ('Filter out SNPs with inconsistent direction:\n-----------------')
    logger.info ('{} SNPs with inconsistent direction'.format(np.sum(Idx==False)))
    logger.info ('\n save invalid SNPs to "{}"'.format(outfile))
    logger.info('\n-------------------------------------------')
    tmpD = sumDat.loc[Idx==False,:]
    tmpD.to_csv(outfile, na_rep='NA', compression='gzip', index=False, sep='\t')
    return (sumDat.loc[Idx==True])

