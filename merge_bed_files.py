# Merge together several .bed files, ignoring all potential warnings from plink.
#
# To run this tool:
#    python merge_bed_files.py --bed ~/1000Genome/phase3/build37_released/*.bed --out merged

import argparse
import glob
import itertools
import os.path
import os
import subprocess
import sys
import pandas as pd

def parse_args(args):
    parser = argparse.ArgumentParser(description="Merge together several .bed files")
    parser.add_argument("--bed", type=str, help="Filename of input .bed file, or pattern (for example '~/1000Genome/phase3/build37_released/*.bed')")
    parser.add_argument("--out", default=r"merged", type=str, help="Filename of output .bed file (without extention")
    return parser.parse_args(args)

def execute_command(command):
    print("Execute command: {}".format(command))
    print(subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].decode("utf-8"))
    #print(subprocess.check_output(command.split()).decode("utf-8"))

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
    execute_command('plink --memory 4096 --bfile {0} --merge-list mergelist.txt --allow-no-sex --make-bed --out {1}'.format(first, output_bfile))
    os.remove('mergelist.txt')

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])

    # Find all .bed filenames (without extention)
    files = [os.path.splitext(file)[0] for file in glob.glob(args.bed)]

    merge(files, args.out)
    missnp_file = '{0}-merge.missnp'.format(args.out)
    if os.path.exists(missnp_file):
        # Handle merge failure as described here: https://www.cog-genomics.org/plink2/data#merge3
        for file in files:
            exclude_snps(file, missnp_file, '{0}.filter'.format(file))
        merge(['{0}.filter'.format(file) for file in files], args.out)
        for file in files: map(os.remove, glob.glob('{0}.filter'.format(file)))

    print("Done.")
