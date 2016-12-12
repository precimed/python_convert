# make_ld_matrix.py

Make LD matrix from reference data. Output either in matlab format or as dense lower triangular text file.
To run this tool you need to download reference data from http://ctg.cncr.nl/software/magma (for example g1000_eur).

For info run `python make_ld_matrix.py --help`.


## Usage

### only abel:
```
qlogin --account=nn9114k --mem-per-cpu=2000  --cpus-per-task=2

module load python2
module load plink2
module load octave

pip install --user pandas
pip install --user scipy

cd /work/users/$USER/
```

### general:

Download data:
  ```
(
  mkdir data10g
  cd data10g
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr[0-9]*vcf.gz
)
  ```
  
Download code:
  ```
git clone https://github.com/precimed/python_convert
cd python_convert/make_ld_matrix
  ```

Set path to your reference file:
  ```
  reffile=/work/users/$USER/2558411_ref.bim 
  ```

Create Matlab matrices:
  ```
for i in ../../data10g/*.vcf.gz; do
  python make_ld_matrix.py --vcf  $i \
    --ref $reffile  --savemat tmp/$(basename $i).map \
    --plink 'plink --memory 3600 --threads 2';
done
  ```

Convert to Matlab sparse matrices:
```
for f in tmp/*.map.mat; do
  # repace ".map.mat" with ".sparse.mat"
  outfile=${f%.map.mat}.sparse.mat
  # matlab script
  mscript="
  load $f
LDmat = sparse(double(id1),double(id2),true,double(nsnp),double(nsnp));
LDmat = LDmat | speye(double(nsnp));
LDmat = LDmat | (LDmat - LDmat');
save(\"$outfile\", 'LDmat', '-v7.3')
"
  # run matlab script
  echo "$mscript" | octave --silent
done
```
Inside the folder "./tmp/" you now have a file "*.sparse.mat" for each chromosome. 
