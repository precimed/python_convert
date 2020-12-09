import pandas as pd
import numpy as np


def get_byte_map():
    """
    Construct mapping between bytes 0..255 and 4-element arrays of a1 genotypes
    from plink bed file.
    Return 256 x 4 array A, where A[i] = [a1, a2, a3, a4], each ai is from {2, -1, 1, 0}.
    """
    genotype_codes = np.array([2, -1, 1, 0],dtype=np.int8)
    byte_map = np.empty((256,4), dtype=np.int8)
    for b in range(256):
        for a in range(4):
            byte_map[b,a] = genotype_codes[(b >> 2*a) & 3]
    return byte_map


class Plink(object):
    bim = None # pd.DataFrame
    fam = None # pd.DataFrame
    bed = None # np.memmap(dtype=np.uint8), this also contain extra bits if n_samples%4 != 0

    byte_map = get_byte_map()

    def __init__(self, bfile):
        print(f"Loading plink {bfile}")
        self._load_bim(bfile)
        self._load_fam(bfile)
        self._load_bed(bfile)
    
    def _load_bim(self, bfile):
        self.bim = pd.read_csv(f'{bfile}.bim', sep='\t', header=None,
                               names=["chr","snp","cm","bp","a1","a2"])
    
    def _load_fam(self, bfile):
        self.fam = pd.read_csv(f'{bfile}.fam',sep='\t',header=None,
                               names=["fid","iid","father_id","mother_id","sex","pheno"])

    def _load_bed(self, bfile):
        if self.bim is None:
            self._load_bim(bfile)
        if self.fam is None:
            self._load_fam(bfile)
        bedf = f'{bfile}.bed'
        magic_bits = np.fromfile(bedf, count=3, dtype=np.uint8)
        if (magic_bits != [108,27,1]).any():
            raise Exception(f"{bedf} file is not a valid bed file!")
        n_snps = self.bim.shape[0]
        n_samples = self.fam.shape[0]
        n_cols = n_samples//4
        if 4*n_cols != n_samples:
            n_cols += 1
        self.bed = np.memmap(bedf, dtype=np.uint8, offset=3, mode='r', shape=(n_snps,n_cols))
        
    def get_geno(self, snp_ii=None):
        """
        Get genotypes for SNPs with indices from snp_ii.
        Args:
            snp_ii  : np.array of SNP indices
        """
        n_snps = self.bim.shape[0]
        n_samples = self.fam.shape[0]
        if snp_ii is None:
            snp_ii = np.arange(n_snps)
        assert max(snp_ii) < n_snps, f"SNP index cannot be > {n_snps-1}"
        n_cols = 4*(n_samples//4)
        if n_cols != n_samples:
            n_cols += 4
        samp_geno = self.byte_map[self.bed[snp_ii]].reshape((len(snp_ii),n_cols))
        return samp_geno[:,:n_samples]


if __name__ == '__main__':
    # Example:
    bfile = '/path/to/plink_bfile'
    plink = Plink(bfile)

    # read genotypes of 0-th and 10-th variants
    geno_arr = plink.get_geno([0,10]) 
    print(geno_arr.shape)

    # read genotypes of all variants in the bfile
    geno_arr = plink.get_geno() 
    print(geno_arr.shape)
