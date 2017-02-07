import sys, os, re, logging, datetime
import numpy as np
import pandas as pd
import gzip
from six import iteritems, string_types
from pkg_resources import parse_version
from itertools import islice

pdlow = parse_version(pd.__version__) < parse_version('0.17.0')

def _get_str_list_sign(str_list):
    return np.array([-1 if e[0]=='-' else 1 for e in str_list], dtype=np.int)

null_values = {

    'LOG_ODDS': 0,
    'BETA': 0,
    'OR': 1,
    'Z': 0
}

default_cnames = {

    # RS NUMBER
    'SNP': 'SNP',
    'MARKERNAME': 'SNP',
    'SNPID': 'SNP',
    'RS': 'SNP',
    'RSID': 'SNP',
    'RS_NUMBER': 'SNP',
    'RS_NUMBERS': 'SNP',
    # CHROMOSOME
    'CHR': 'CHR',
    'CHROMOSOME' : 'CHR',
    # POSITION
    'POS': 'POS',
    'BP': 'POS',
    'POSITION' : 'POS',
    # NUMBER OF STUDIES
    'NSTUDY': 'NSTUDY',
    'N_STUDY': 'NSTUDY',
    'NSTUDIES': 'NSTUDY',
    'N_STUDIES': 'NSTUDY',
    # P-VALUE
    'P': 'P',
    'PVALUE': 'P',
    'P_VALUE':  'P',
    'PVAL': 'P',
    'P_VAL': 'P',
    'GC_PVALUE': 'P',
    # ALLELE 1
    'A1': 'A1',
    'ALLELE1': 'A1',
    'ALLELE_1': 'A1',
    'EFFECT_ALLELE': 'A1',
    'REFERENCE_ALLELE': 'A1',
    'INC_ALLELE': 'A1',
    'EA': 'A1',
    # ALLELE 2
    'A2': 'A2',
    'ALLELE2': 'A2',
    'ALLELE_2': 'A2',
    'OTHER_ALLELE': 'A2',
    'NON_EFFECT_ALLELE': 'A2',
    'NON_EFF_ALLELE': 'A2',
    'DEC_ALLELE': 'A2',
    'NEA': 'A2',
    # N
    'N': 'N',
    'NCASE': 'N_CAS',
    'CASES_N': 'N_CAS',
    'N_CASE': 'N_CAS',
    'N_CASES': 'N_CAS',
    'N_CONTROLS': 'N_CON',
    'N_CAS': 'N_CAS',
    'N_CON': 'N_CON',
    'N_CASE': 'N_CAS',
    'NCONTROL': 'N_CON',
    'CONTROLS_N': 'N_CON',
    'N_CONTROL': 'N_CON',
    'WEIGHT': 'N',  # metal does this. possibly risky.
    # SIGNED STATISTICS
    'ZSCORE': 'Z',
    'Z-SCORE': 'Z',
    'GC_ZSCORE': 'Z',
    'Z': 'Z',
    'OR': 'OR',
    'B': 'BETA',
    'BETA': 'BETA',
    'LOG_ODDS': 'LOG_ODDS',
    'EFFECTS': 'BETA',
    'EFFECT': 'BETA',
    'SIGNED_SUMSTAT': 'SIGNED_SUMSTAT',
    # STANDARD ERROR
    'SE' : 'SE',
    # INFO
    'INFO': 'INFO',
    # MAF
    'EAF': 'FRQ',
    'FRQ': 'FRQ',
    'MAF': 'FRQ',
    'FRQ_U': 'FRQ',
    'F_U': 'FRQ',
}

describe_cname = {
    'SNP': 'Variant ID (e.g., rs number)',
    'CHR': 'Chromosome number',
    'POS': 'Base-pair position',
    'P': 'p-Value',
    'A1': 'Allele 1, interpreted as ref allele for signed sumstat.',
    'A2': 'Allele 2, interpreted as non-ref allele for signed sumstat.',
    'N': 'Sample size',
    'N_CAS': 'Number of cases',
    'N_CON': 'Number of controls',
    'Z': 'Z-score (0 --> no effect; above 0 --> A1 is trait/risk increasing)',
    'OR': 'Odds ratio (1 --> no effect; above 1 --> A1 is risk increasing)',
    'BETA': '[linear/logistic] regression coefficient (0 --> no effect; above 0 --> A1 is trait/risk increasing)',
    'LOG_ODDS': 'Log odds ratio (0 --> no effect; above 0 --> A1 is risk increasing)',
    'SE': 'standard error of the effect size',
    'INFO': 'INFO score (imputation quality; higher --> better imputation)',
    'FRQ': 'Allele frequency',
    'SIGNED_SUMSTAT': 'Directional summary statistic as specified by --signed-sumstats.',
    'NSTUDY': 'Number of studies in which the SNP was genotyped.',
    'UNKNOWN': 'Unknown column type (will be skipped).'
}

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

def print_header(fh, lines=5):
    (openfunc, compression) = get_compression(fh)
    with openfunc(fh) as f:
        for line in islice(f, lines):
            line = line if isinstance(line, string_types) else line.decode('utf-8') 
            print(line.rstrip('\n'))

def read_header(fh):
    '''Read the first line of a file and returns a list with the column names.'''
    (openfunc, compression) = get_compression(fh)
    firstline = openfunc(fh).readline()
    if not isinstance(firstline, string_types):
        firstline = firstline.decode('utf-8')
    return [x.rstrip('\n') for x in firstline.split()]

def get_cname_map(flag, default, ignore):
    '''
    Figure out which column names to use.
    Priority is
    (1) ignore everything in ignore
    (2) use everything in flags that is not in ignore
    (3) use everything in default that is not in ignore or in flags
    The keys of flag are cleaned. The entries of ignore are not cleaned. The keys of defualt
    are cleaned. But all equality is modulo clean_header().
    '''
    clean_ignore = [clean_header(x) for x in ignore]
    cname_map = {x: flag[x] for x in flag if x not in clean_ignore}
    cname_map.update(
        {x: default[x] for x in default if x not in clean_ignore + list(flag)})
    return cname_map

def get_compression(fh):
    '''
    Read filename suffixes and figure out whether it is gzipped,bzip2'ed or not compressed
    '''
    if fh.endswith('gz'):
        compression = 'gzip'
        openfunc = gzip.open
    elif fh.endswith('bz2'):
        compression = 'bz2'
        openfunc = bz2.BZ2File
    else:
        openfunc = open
        compression = None

    return openfunc, compression

def clean_header(header):
    '''
    For cleaning file headers.
    - convert to uppercase
    - replace dashes '-' with underscores '_'
    - replace dots '.' (as in R) with underscores '_'
    - remove newlines ('\n')
    '''
    return header.upper().replace('-', '_').replace('.', '_').replace('\n', '')

def parse_flag_cnames(cname_options, signed_sumstats):
    '''
    Parse flags that specify how to interpret nonstandard column names.
    flag_cnames is a dict that maps (cleaned) arguments to internal column names
    '''
    flag_cnames = {clean_header(x[0]): x[1]
                   for x in cname_options if x[0] is not None}
    null_value = None
    if signed_sumstats is not None:
        try:
            cname, null_value = signed_sumstats.split(',')
            null_value = float(null_value)
            flag_cnames[clean_header(cname)] = 'SIGNED_SUMSTAT'
        except ValueError:
            #log.log('The argument to --signed-sumstats should be column header comma number.')
            raise

    return [flag_cnames, null_value]

def make_metal_script(files, out='meta'):
    with open('metal_script.txt', 'w') as f:
        f.write('SCHEME STDERR\n')
        f.write('SEPARATOR WHITESPACE\n')

        for file in files:
            cname_translation = find_column_name_translation(sumstats=file)
            cname_description = {x: describe_cname[cname_translation[x]] for x in cname_translation if cname_translation[x] != 'UNKNOWN'}
            print('Interpreting column names from {0} as follows:'.format(file))
            print('\n'.join([x + ':\t' + cname_description[x] for x in cname_description]) + '\n')
            cname_skip = [x for x in cname_translation if cname_translation[x ] == 'UNKNOWN']
            if cname_skip: print('Skip the remaining columns ({}).'.format(', '.join(cname_skip)))
            cname = {v: k for k, v in iteritems(cname_translation)}

            se = cname.get('SE', 'UNKNOWN_COLUMN')
            pvalue = cname.get('P', 'UNKNOWN_COLUMN')
            snp = cname.get('SNP', 'UNKNOWN_COLUMN')
            A1 = cname.get('A1', 'UNKNOWN_COLUMN')
            A2 = cname.get('A2', 'UNKNOWN_COLUMN')
            effect = 'UNKNOWN_COLUMN'
            if 'BETA' in cname: effect = cname['BETA']
            elif 'LOG_ODDS' in cname: effect = cname['LOG_ODDS']
            elif 'OR' in cname: effect = 'log({0})'.format(cname['OR'])
            if effect == 'UNKNOWN_COLUMN': print('WARNING: Effect size column not detected in {0}'.format(file))
            if snp == 'UNKNOWN_COLUMN': print('WARNING: SNP column not detected in {0}'.format(file))
            if A1 == 'UNKNOWN_COLUMN' or A2 == 'UNKNOWN_COLUMN': print('WARNING: A1/A2 columns not detected in {0}'.format(file))
            if se == 'UNKNOWN_COLUMN': print('WARNING: SE column not detected in {0}'.format(file))
            if pvalue == 'UNKNOWN_COLUMN': print('WARNING: P value column not detected in {0}'.format(file))

            f.write('MARKER {0}\n'.format(snp))
            f.write('ALLELE {0} {1}\n'.format(A1, A2))
            f.write('EFFECT {0}\n'.format(effect))
            f.write('STDERR {0}\n'.format(se))
            f.write('PVALUE {0}\n'.format(pvalue))
            f.write('PROCESS {0}\n'.format(file))

        f.write('OUTFILE {0} .tbl\n'.format(out))
        f.write('ANALYZE\n')
        f.write('QUIT')
    print('Done. Now you should run "metal script metal_script.txt"')

def find_column_name_translation(sumstats, snp=None, chromosome=None, pos=None,
                                 a1=None, a2=None, a1_inc=None, p=None,
                                 frq=None, info=None,
                                 signed_sumstats=None, odds_ratio=None, beta=None,
                                 ignore=None, daner=False,
                                 nstudy=None, N_col=None, N_cas_col=None, N_con_col=None):
    '''
    Detects naming of column in summary stats file.
    Returns a dictionary where key corresponds to column name from summary stat file, and
    value is one of the following:
    SNP, CHR, POS, P, A1, A2, N, N_CAS, N_CON, Z, OR, BETA, LOG_ODDS, INFO, FRQ, SIGNED_SUMSTAT, or NSTUDY.
    Corresponding parameters of this function allow to provide a custom name that could not be detected automaticaly.
    The code is largely copied from https://github.com/bulik/ldsc/blob/master/munge_sumstats.py
    '''
    cname_options = [
        [nstudy, 'NSTUDY', '--nstudy'],
        [snp, 'SNP', '--snp'],
        [chromosome, 'CHR', '--chr'],
        [odds_ratio, 'OR', '--or'],
        [beta, 'BETA', '--beta'],
        [pos, 'POS', '--pos'],
        [N_col, 'N', '--N'],
        [N_cas_col, 'N_CAS', '--N-cas-col'],
        [N_con_col, 'N_CON', '--N-con-col'],
        [a1, 'A1', '--a1'],
        [a2, 'A2', '--a2'],
        [p, 'P', '--P'],
        [frq, 'FRQ', '--nstudy'],
        [info, 'INFO', '--info']
    ]

    file_cnames = read_header(sumstats)  # note keys not cleaned
    flag_cnames, signed_sumstat_null = parse_flag_cnames(cname_options, signed_sumstats)
    if ignore:
        ignore_cnames = [clean_header(x) for x in ignore.split(',')]
    else:
        ignore_cnames = []

    # remove LOG_ODDS, BETA, Z, OR from the default list
    if signed_sumstats is not None or a1_inc:
        mod_default_cnames = {x: default_cnames[x] for x in default_cnames if default_cnames[x] not in null_values}
    else:
        mod_default_cnames = default_cnames

    cname_map = get_cname_map(flag_cnames, mod_default_cnames, ignore_cnames)

    if daner:
        frq_u = next(iter(filter(lambda x: x.startswith('FRQ_U_'), file_cnames)))
        frq_a = next(iter(filter(lambda x: x.startswith('FRQ_A_'), file_cnames)))
        N_cas = float(frq_a[6:])
        N_con = float(frq_u[6:])
        #log.log('Inferred that N_cas = {N1}, N_con = {N2} from the FRQ_[A/U] columns.'.format(N1=N_cas, N2=N_con))

        # drop any N, N_cas, N_con or FRQ columns
        for c in ['N', 'N_CAS', 'N_CON', 'FRQ']:
            for d in [x for x in cname_map if cname_map[x] == 'c']:
                del cname_map[d]

        cname_map[frq_u] = 'FRQ'

    cname_translation = {x: cname_map.get(clean_header(x), 'UNKNOWN') for x in file_cnames}  # note keys not cleaned

    # In the original version the code here performed the following validation
    # 1. One and only one signed sumstat column is found (unless a1_inc is not None)
    # 2. 'SNP' and 'P' columns are present in the data
    # 3. Alleles 'A1' and 'A2' are present

    return cname_translation
