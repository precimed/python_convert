# Download reference data from http://ctg.cncr.nl/software/magma (for example g1000_eur)
# Then you can run the tool as follows:
#    python make_ld_matrix.py --ref 2558411_ref.bim --bfile g1000_eur --ld_window_r2 0.1 --savemat ldmat_p1.mat


from subprocess import call, check_output
import subprocess
import pandas as pd
import argparse
import sys

def execute_command(command):
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(process.communicate()[0].decode("utf-8"))
    #print(subprocess.check_output(command.split()).decode("utf-8"))


def parse_args(args):
    parser = argparse.ArgumentParser(description="Generate LD matrix from genotype matrix")
    parser.add_argument("--ref", type=str, help="Reference file (for example 2558411_ref.bim or 9279485_ref.bim.")
    parser.add_argument("--bfile", type=str, help="Genotypes in plink binary format")
    parser.add_argument("--ld_window_kb", default=5000, type=int, help="Window in KB")
    parser.add_argument("--ld_window_r2", default=0.1, type=float, help="LD r2 threshold")
    parser.add_argument("--chunksize", default=100000, type=int, help="Chunk size when reading ld matrix")
    parser.add_argument("--plink", default="plink", type=str, help="location of plink executable")
    parser.add_argument("--savemat", default=None, type=str, help="Generate matfile for Matlab.")
    parser.add_argument("--saveltm", default=None, type=str, help="Generate 'ltm' --- lower triangular matrix in plain text format.")
    return parser.parse_args(args)


def make_ld_matrix(args):
    if not args.savemat and not args.saveltm:
        raise ValueError('No output requested, use --savemat or --saveltm')

    # Read the template
    ref = pd.read_csv(args.ref, delim_whitespace=True)
    nsnp = ref.shape[0]
    rs_to_id = dict(zip(ref['SNP'], ref.index))
    if len(rs_to_id) != nsnp: raise ValueError("Duplicated SNPs found in the reference file")

    # Filter according to the template, then create LD file in table format
    execute_command('{0} --bfile {1} --make-bed --out tmp --extract {2}'.format(args.plink, args.bfile, args.ref))
    execute_command('{0} --bfile tmp --r2 --ld-window-kb {1} --ld-window 999999 --ld-window-r2 {2} --out tmp'.format(args.plink, args.ld_window_kb, args.ld_window_r2))

    # Read resulting LD matrix
    reader = pd.read_csv('tmp.ld', delim_whitespace=True, chunksize=args.chunksize)

    id1 = []; id2 = []; val = []
    for i, df in enumerate(reader):
        id1.extend(list(df['SNP_A'].map(lambda x: rs_to_id[x])))
        id2.extend(list(df['SNP_B'].map(lambda x: rs_to_id[x])))
        val.extend(list(df['R2']))

    # Output the result as lower diagonal matrix
    if args.saveltm:
        from scipy.sparse import csr_matrix
        assert(all([(i < j) for (i, j) in zip(id1, id2)]))  # expect that plink output lower diagonal matrix
        csr = csr_matrix((val, (id2, id1)), shape=(nsnp, nsnp))

        with open(args.saveltm, 'w') as result:
            result.write('1.0\n')
            for i in range(1, nsnp):
                values =  csr[i, :].todense()[0, 0:i].A1
                values_str = '\t'.join(str(x) for x in values)
                result.write('{0}\t1.0\n'.format(values_str))
    
    # Output the result in matlab format
    if args.savemat:
        import scipy.io as sio
        sio.savemat(
            args.savemat, {'id1':id1, 'id2':id2, 'nsnp':nsnp},
            format='5', do_compression=False, oned_as='column')

        print("""
The results are saved into {0}. Now you should open matlab and execute the following commands to re-save the result as matlab sparse matrix:
    load {0}
    LDmat = sparse(double(1+id1),double(1+id2),true,double(nsnp),double(nsnp));
    LDmat = LDmat | speye(double(nsnp));
    LDmat = LDmat | (LDmat - LDmat');
    save('LDmat.mat', 'LDmat', '-v7.3')
""".format(args.savemat))


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    make_ld_matrix(args)
    print("Done.")
