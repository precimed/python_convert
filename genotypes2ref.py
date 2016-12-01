# Align genotypes to reference file, e.g.
# - Extract SNPs from the reference (merge by CHR:POS --- require data to be on the same genomic build)
# - Extract subset of individuals (for example, european population)
# - Merge SNPs together (for example if input data is split by chromosome)
#
# To run this tool:
# - Download *.vcf.gz files from 1000 Genome project ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
# - Run the tool as follows 
#    python genotypes2ref.py --dir 1000Genome/phase3/build37_released --ref 2558411_ref.bim

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
    parser.add_argument("--dir", type=str, help="Folder with input *.vcf.gz files")
    parser.add_argument("--keep", default=r"data/EUR_subj.list", type=str, help="Extract SNPs and keep only EUR individuals")
    return parser.parse_args(args)

def execute_command(command):
    print("Execute command: {}".format(command))
    print(subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].decode("utf-8"))
    #print(subprocess.check_output(command.split()).decode("utf-8"))

def process_vcf_file(vcf_file, df_ref, keep_file):
    [_, filename] = os.path.split(vcf_file)
    output_dir = os.path.join(working_path, 'tmp')
    bfile = os.path.join(working_path, 'tmp', filename)
    snpidlist = os.path.join(working_path, 'tmp', filename + '.snpidlist.txt')
    join_file = os.path.join(working_path, 'tmp', filename + '.joined')

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

def exclude_snps(bfile_in, snps_file, bfile_out):
    execute_command('plink --memory 4096 --bfile {0} --exclude {1} --make-bed --out {2}'.format(bfile_in, snps_file, bfile_out))
    
def merge(files, output_bfile):
    missnp_file = '{0}-merge.missnp'.format(output_bfile)
    if os.path.exists(missnp_file):
        os.remove(missnp_file)
    first = files[0]
    with open('mergelist.txt', 'w') as mergelist:
        for filename in files[1:]:
            mergelist.write('{0}.bed {0}.bim {0}.fam\n'.format(filename))
    execute_command('plink --memory 4096 --bfile {0} --merge-list mergelist.txt --allow-no-sex --freq --make-bed --out {1}'.format(first, output_bfile))
    os.remove('mergelist.txt')

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])

    vcf_files=[file for file in glob.glob(os.path.join(working_path, '*.vcf.gz'))
               if ('chrX' not in file) and ('chrY' not in file)]
    assert len(vcf_files) == 22
        
    # Read reference file
    df_ref = pd.read_csv(args.ref, delim_whitespace=True)
    assert df_ref.duplicated(['CHR', 'BP']).sum() == 0
    assert df_ref.duplicated(['SNP']).sum() == 0

    for vcf_file in vcf_files: process_vcf_file(vcf_file, df_ref, args.keep)

    # Find all .bed filenames (without extention)
    files = [os.path.splitext(file)[0] for file in glob.glob(os.path.join(args.dir, 'tmp', '*.joined.bed'))]

    output_bfile = os.path.join(args.dir, 'tmp', 'merged')
    merge(files, output_bfile)
    missnp_file = '{0}-merge.missnp'.format(output_bfile)
    if os.path.exists(missnp_file):
        # Handle merge failure as described here: https://www.cog-genomics.org/plink2/data#merge3
        for file in files:
            exclude_snps(file, missnp_file, '{0}.filter'.format(file))
        merge(['{0}.filter'.format(file) for file in files], output_bfile)

    print("Done.")
