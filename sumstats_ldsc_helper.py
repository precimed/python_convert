# Examples:
# python sumstats_ldsc_helper.py annot PGC_SCZ_2014.csv PGC_SCZ_2014_ldscores/{}.annot.gz --annot 1000G_Phase3_baseline_ldscores/baseline.{}.annot.gz  

import pandas as pd
import numpy as np
import os.path
import sys
import argparse
from intervaltree import IntervalTree

def parse_args(args):
    parser = argparse.ArgumentParser(description="Miscellaneous utilities to work with statistics and LDSC regression method")
    subparsers = parser.add_subparsers()
    parser_annot = subparsers.add_parser("annot", help="Create binary annotations from 0.1, 0.01 and 0.001 p-value stratuums of the summary statistic")
    parser_annot.add_argument("sumstats_file", type=str, help="Input file with summary statistics")
    parser_annot.add_argument("output_file", type=str, help="Path to the output file to place the results")
    parser_annot.add_argument("--annot", type=str, help="Path to baseline.CHR.annot.gz file from 1000G_Phase3_baseline_ldscores")
    parser_annot.add_argument("--force", action="store_true", default=False, help="Force overwrite target files if they exist.")
    parser_annot.add_argument("--window", default=2, type=int, help="Window to include into the binary annotation around each SNP that pass p-value threshold.")
    parser_annot.set_defaults(func=make_annot)
    return parser.parse_args(args)

### =================================================================================
###                          Implementation for parser_annot
### ================================================================================= 
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

def make_annot(args):
    """
    Create binary annotations from 0.1, 0.01 and 0.001 p-value stratuums of the summary statistic.
    """
    check_input_file(args.sumstats_file)
    for chri in range(1, 23): check_output_file(args.output_file.format(chri), args.force)

    print('Reading summary statistics file {}...'.format(args.sumstats_file))
    sumstats = pd.read_csv(args.sumstats_file, delim_whitespace=True, usecols=['PVAL', 'CHR', 'BP'])
    print('Done, read {} SNPs.'.format(sumstats.shape[0]))

    for chri in range(1, 23):
        print('Processing chromosome {}...'.format(chri))
        df = pd.read_csv(args.annot.format(chri), delim_whitespace=True)
        df = df[['CHR', 'BP', 'SNP', 'CM']].copy()
        for pthresh, label in [(0.1, '.1'), (0.01, '.01'), (0.001, '.001')]:
            sumstatsCHR = sumstats[sumstats.CHR == chri].copy(deep=True)
            print('{} markers, {} of them are on chr {}, {} of them have p-value below {}'.format(sumstats.shape[0], sumstatsCHR.shape[0], chri, (sumstatsCHR.PVAL < pthresh).sum(), pthresh))
            itree = IntervalTree.from_tuples(zip(sumstatsCHR[sumstatsCHR.PVAL < pthresh].BP - args.window, sumstatsCHR[sumstatsCHR.PVAL < pthresh].BP + args.window))
            itree.merge_overlaps()
            print('Found {} intervals, average length {}'.format(len(itree),  sum([i.length() for i in itree])/len(itree)))

            annot_binary = [int(bool(itree[p])) for p in df.BP]
            df['PVAL{}'.format(label)] = annot_binary
            print('{} markers out of {} ({}%) belongs to the annotation'.format(sum(annot_binary), len(annot_binary), 100 * sum(annot_binary) / len(annot_binary)))
        df.to_csv(args.output_file.format(chri), index=False, sep='\t', compression='gzip')
        print('Results saved to {}'.format(args.output_file.format(chri)))

### =================================================================================
###                                Main section
### ================================================================================= 
if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    args.func(args)
    print("Done")
