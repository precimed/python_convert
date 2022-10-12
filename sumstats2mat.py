#!/usr/bin/env python
import os
import sys
import pandas as pd
import argparse
import scipy.io as sio
import numpy as np
from collections import namedtuple

Cols = namedtuple('Cols', ['SNP',  'PVAL', 'A1',           'A2',          'N', 'NCASE', 'NCONTROL', 'Z'])
cols = Cols._make(        ['RSID', 'P',    'EffectAllele', 'OtherAllele', 'N', 'CaseN', 'ControlN', 'Z'])

__version__ = '1.0.0'
MASTHEAD = "***********************************************************************\n"
MASTHEAD += "* sumstats2mat.py: Converts summary statistics from csv to matlab format\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C) Norwegian Centre for Mental Disorders Research / University of Oslo\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "***********************************************************************\n"

_base_complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
def _reverse_complement_variant(variant):
    # variant should be a 2-elemet sequence with upper case string elements
    return ("".join([_base_complement[b] for b in variant[0][::-1]]),
            "".join([_base_complement[b] for b in variant[1][::-1]]))

def check_input_file(file):
    if not os.path.isfile(file):
        raise ValueError("Input file does not exist: {f}".format(f=file))

def check_output_file(file, force=False):
    # Delete target file if user specifies --force option
    if force:
        try:
            os.remove(file)
        except OSError:
            pass

    # Otherwise raise an error if target file already exists
    if os.path.isfile(file) and not force:
        raise ValueError("Output file already exists: {f}".format(f=file))

    # Create target folder if it doesn't exist
    output_dir = os.path.dirname(file)
    if output_dir and not os.path.isdir(output_dir): os.makedirs(output_dir)  # ensure that output folder exists

def make_mat(args):
    if args.out is None: raise ValueError('--out is required.')

    check_input_file(args.ref)
    check_input_file(args.sumstats)
    check_output_file(args.out, args.force)

    reader = pd.read_csv(args.sumstats, delim_whitespace=True, chunksize=args.chunksize, float_precision='high')
    df_out = None
    for chunk_index, ss_chunk in enumerate(reader):
        # (BEGIN) special handling of the first chunk
        if chunk_index==0:
            columns = list(ss_chunk.columns)

            required_cols = [cols.SNP, cols.A1, cols.A2, cols.Z]
            if (set(required_cols) - set(columns)):
                absent_cols = set(required_cols) - set(columns)
                err_msg = ("Columns {} are missing from the --sumstats file {}").format(', '.join(absent_cols), args.sumstats)
                raise(RuntimeError(err_msg))

            n_col = cols.N if cols.N in columns else None
            ncase_col = cols.NCASE if cols.NCASE in columns else None
            ncontrol_col = cols.NCONTROL if cols.NCONTROL in columns else None
            if (not args.without_n) and ((n_col is None) and ((ncase_col is None) or (ncontrol_col is None))):
                raise(ValueError('Sample size column is not detected in {}. Expact either N or NCASE, NCONTROL column.'.format(args.sumstats)))

            print('Reading reference file {}...'.format(args.ref))
            ref_reader = pd.read_csv(args.ref, sep='\t', usecols=['SNP', 'A1', 'A2'], chunksize=args.chunksize)
            ref_dict = {}
            for ref_chunk in ref_reader:
                ref_chunk.drop(ref_chunk.index[np.logical_not(ref_chunk['A1'].str.upper().str.match('^[ACTG]*$')) | np.logical_not(ref_chunk['A2'].str.upper().str.match('^[ACTG]*$'))], inplace=True)
                if ref_chunk.empty: continue
                gtypes = zip(ref_chunk['A1'].apply(str.upper),ref_chunk['A2'].apply(str.upper))
                #TODO?: add check whether some id is already in ref_dict
                ref_dict.update(dict(zip(ref_chunk['SNP'], gtypes)))
            ref_dict = {i: (variant, _reverse_complement_variant(variant),
                            variant[::-1], _reverse_complement_variant(variant[::-1]))
                        for i, variant in ref_dict.items()}
            ref_snps = pd.read_csv(args.ref, sep='\t', usecols=['SNP'], squeeze=True)
            #TODO?: add check whether ref_snps contains duplicates
            print("Reference dict contains {d} snps.".format(d=len(ref_dict)))

            print('Reading summary statistics file {}...'.format(args.sumstats))
            print('Column types: ' + ', '.join([column + ':' + str(dtype) for (column, dtype) in zip(ss_chunk.columns, ss_chunk.dtypes)]))
        # (END) special handling of the first chunk

        ss_chunk = ss_chunk.loc[ss_chunk[cols.SNP].isin(ref_dict),:]
        if ss_chunk.empty: continue
        gtypes = list(zip(ss_chunk[cols.A1].apply(str.upper),ss_chunk[cols.A2].apply(str.upper)))
        # index of SNPs that have the same alleles as indicated in reference
        ind = [gt in ref_dict[sid] for sid, gt in zip(ss_chunk[cols.SNP], gtypes)]
        ss_chunk = ss_chunk.loc[ind,:]
        gtypes = [gt for gt, j in zip(gtypes, ind) if j]
        log10pv = -np.log10(ss_chunk[cols.PVAL].values)
        # not_ref_effect = [
        #   1 if effect allele in data == other allele in reference
        #   -1 if effect allele in data == effect allele in reference ]
        # So zscores with positive effects will be positive and zscores with
        # negative effects will stay negative, since
        # stats.norm.ppf(ss_chunk[cols.PVAL]*0.5) is always negetive (see zvect
        # calculation below).
        not_ref_effect = np.array([1 if gt in ref_dict[sid][:2] else -1
            for sid, gt in zip(ss_chunk[cols.SNP], gtypes)])
        #TODO: check proportion of positive and negative effects
        zvect = ss_chunk[cols.Z].values*not_ref_effect
        ind_ambiguous = [j for j,gt in enumerate(gtypes) if gt == _reverse_complement_variant(gt)[::-1]]
        # set zscore of ambiguous SNPs to nan
        zvect[ind_ambiguous] = np.nan
        #TODO: check whether output df contains duplicated rs-ids (warn)

        # reindex by SNP, add required columns and drop unnecessary columns
        ss_chunk.index = ss_chunk[cols.SNP]
        # add required columns
        ss_chunk["logpvec"] = log10pv
        ss_chunk["zvec"] = zvect
        if not args.without_n:
            if n_col is None:
                nvec = 4./(1./ss_chunk[ncase_col] + 1./ss_chunk[ncontrol_col])
            else:
                nvec = ss_chunk[n_col].values
            ss_chunk["nvec"] = nvec

        cols2drop = [c for c in ss_chunk.columns if (c not in ['logpvec', 'zvec', 'nvec'])]
        ss_chunk.drop(cols2drop, axis=1, inplace=True)

        if df_out is None:
            df_out = ss_chunk.copy()
        else:
            df_out = df_out.append(ss_chunk)

        print("{f}: {n} lines processed, {m} SNPs matched with reference file".format(f=args.sumstats, n=(chunk_index+1)*args.chunksize, m=len(df_out)))

    if df_out.empty: raise(ValueError("No SNPs match after joining with reference data"))
    dup_index = df_out.index.duplicated(keep=False)
    if dup_index.any():
        print("Duplicated SNP ids detected:")
        print(df_out[dup_index])
        print("Keeping only the first occurance.")
    df_out = df_out[~df_out.index.duplicated(keep='first')]
    # allign index accordind order of SNPs in ref, insert NaN rows for SNPs that
    # present in ref but absent in sumstats file
    df_out = df_out.reindex(ref_snps)

    print('Writing .mat file...')
    save_dict = {c+args.trait: df_out[c].astype(np.float64).values for c in df_out.columns}
    sio.savemat(args.out, save_dict, format='5', do_compression=False,
        oned_as='column', appendmat=False)
    print("%s created" % args.out)

### =================================================================================
###                                Main section
### ================================================================================= 
if __name__ == "__main__":
    parser_mat = argparse.ArgumentParser(description="Create mat files that can "
        "be used as an input for pleiofdr analysis (https://github.com/precimed/pleiofdr/). "
        "Takes a .csv file with summary statistics file as input. The file can be compressed with gzip. "
        "Require columns: SNP, P, A1, A2, Z, and sample size column (either N, or CaseN and ControlN). "
        "Creates corresponding mat files which can be used as an input for pleiofdr analysis. "
        "Only SNPs from the reference file are considered. "
        "Zscores of strand ambiguous SNPs are set to NA. ")

    parser_mat.add_argument("--sumstats", type=str, help="Input file with summary statistics. ")
    parser_mat.add_argument("--ref", type=str, help="[required] Tab-separated file with list of referense SNPs.")
    parser_mat.add_argument("--out", type=str, help="[required] File to output the result. File should end with .mat extension.")
    parser_mat.add_argument("--force", action="store_true", default=False, help="Allow sumstats.py to overwrite output file if it exists.")

    parser_mat.add_argument("--trait", type=str, default='',
        help="Trait name that will be used in mat file. Can be kept empty, in this case the variables will be named 'logpvec', 'zvec' and 'nvec'")
    parser_mat.add_argument("--without-n", action="store_true", default=False,
        help="Proceed without sample size (N or NCASE/NCONTROL)")
    parser_mat.add_argument("--chunksize", default=100000, type=int,
        help="Size of chunk to read the file.")

    args = parser_mat.parse_args()
    make_mat(args)
    
