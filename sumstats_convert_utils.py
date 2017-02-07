import sys, os, re, logging, datetime
import numpy as np
import pandas as pd
import gzip
import six
import itertools

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

def print_header(fh, lines=5):
    (openfunc, compression) = get_compression(fh)
    with openfunc(fh) as f:
        for line in itertools.islice(f, lines):
            line = line if isinstance(line, six.string_types) else line.decode('utf-8') 
            print(line.rstrip('\n'))

def read_header(fh):
    '''Read the first line of a file and returns a list with the column names.'''
    (openfunc, compression) = get_compression(fh)
    firstline = openfunc(fh).readline()
    if not isinstance(firstline, six.string_types):
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
