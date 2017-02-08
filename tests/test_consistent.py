import subprocess
import scipy.io as sio
import numpy as np
import shutil
import os.path

def execute_command(command):
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(process.communicate()[0].decode("utf-8"))
    #print(subprocess.check_output(command.split()).decode("utf-8"))

def run(filename, matfile, effect):
    reffile = r'tests/1234_ref.bim'
    if os.path.isdir('TEMP_FOLDER'): shutil.rmtree('TEMP_FOLDER')
    execute_command(r'python sumstats_convert.py csv {} TEMP_FOLDER/TEST.csv --auto --force'.format(filename))
    execute_command(r'python sumstats_convert.py mat {} TEMP_FOLDER/TEST.csv --traits test --force --effect {}'.format(reffile, effect))

    f1 = sio.loadmat(matfile)
    f2 = sio.loadmat('TEMP_FOLDER/TEST.mat')
    assert(all(np.isfinite(f1['logpvec_test']) == np.isfinite(f2['logpvec_test'])))
    assert(all(np.isfinite(f1['zvec_test']) == np.isfinite(f2['zvec_test'])))
    assert(max(abs(f1['logpvec_test'] - f2['logpvec_test'])) < 1e-10)
    assert(max(abs(f1['zvec_test'] - f2['zvec_test'])) < 1e-10)
    shutil.rmtree('TEMP_FOLDER')

def test01(): run('tests/case01.txt', 'tests/case01.mat', effect='BETA')
def test01gz(): run('tests/case01.txt.gz', 'tests/case01.mat', effect='BETA')
#def test02(): run('tests/case02.txt')
#def test03(): run('tests/case03.txt')
def test04(): run('tests/case04.txt', 'tests/case04.mat', effect='OR')
def test04gz(): run('tests/case04.txt.gz', 'tests/case04.mat', effect='OR')
