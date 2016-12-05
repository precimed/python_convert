# Align genotypes to reference file, e.g.
# - Extract SNPs from the reference (merge by CHR:POS --- require data to be on the same genomic build)
# - Extract subset of individuals (for example, european population)
# - Merge SNPs together (for example if input data is split by chromosome)
#
# To run this tool:
# - Download *.vcf.gz files from 1000 Genome project ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
# - Run the tool as follows 
#    python genotypes2ref.py --vcf ~/1000Genome/phase3/build37_released/*.vcf.gz --ref 2558411_ref.bim

import argparse
import glob
import itertools
import os.path
import os
import subprocess
import sys
import pandas as pd

def parse_args(args):
    parser = argparse.ArgumentParser(description="Generate LD matrix from genotype matrix")
    parser.add_argument("--ref", type=str, help="Reference file (for example 2558411_ref.bim or 9279485_ref.bim.")
    parser.add_argument("--vcf", type=str, help="Filename of input .vcf file, or pattern (for example '~/1000Genome/phase3/build37_released/*.vcf.gz')")
    parser.add_argument("--keep", default=r"data/EUR_subj.list", type=str, help="Extract SNPs and keep only EUR individuals")
    parser.add_argument("--out", default=r"tmp", type=str, help="Folder to output the result")
    return parser.parse_args(args)

def execute_command(command):
    print("Execute command: {}".format(command))
    print(subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].decode("utf-8"))
    #print(subprocess.check_output(command.split()).decode("utf-8"))

def process_vcf_file(vcf_file, df_ref, keep_file, output_dir):
    [_, filename] = os.path.split(vcf_file)
    bfile = os.path.join(output_dir, filename)
    snpidlist = os.path.join(output_dir, filename + '.snpidlist.txt')
    join_file = os.path.join(output_dir, filename)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Convert vcf into bed
    execute_command(r'plink --memory 4096 --vcf {0} --make-bed --out {1}'.format(vcf_file, bfile))

    # Read bim file
    df_bim = pd.read_csv('{}.bim'.format(bfile), header=None, delim_whitespace=True)
    df_bim.columns=['CHR','SNP','GP','POS','A1','A2']

    # Left-merge with reference file by CHR:POS, then output merged RS numberes into snpidlist.txt
    df_bim2 = pd.merge(df_bim, df_ref, how='left', left_on = ['CHR', 'POS'], right_on = ['CHR', 'BP'])
    df_bim2[df_bim2['SNP_y'].notnull()]['SNP_x'].to_csv(snpidlist, index=False)

    # Extract SNPs and keep only EUR individuals
    execute_command(r'plink --memory 4096 --bfile {0} --extract {1} --keep {2} --make-bed --out {3}'.format(bfile, snpidlist, keep_file, join_file))
    return join_file

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])

    vcf_files=[file for file in glob.glob(args.vcf)
               if ('chrX' not in file) and ('chrY' not in file)]
    print(vcf_files)
    # Read reference file
    df_ref = pd.read_csv(args.ref, delim_whitespace=True)
    assert df_ref.duplicated(['CHR', 'BP']).sum() == 0
    assert df_ref.duplicated(['SNP']).sum() == 0

    for vcf_file in vcf_files: process_vcf_file(vcf_file, df_ref, args.keep, args.out)

    print("Done.")
