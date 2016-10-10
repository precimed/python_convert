from __future__ import division
import numpy as np
import scipy.stats as stats
import scipy.special as special
import scipy.linalg as linalg
import scipy.interpolate as interp

# Decorrelation in Statistics: The Mahalanobis Transformation
# Added material to Data Compression: The Complete Reference
# http://www.davidsalomon.name/DC2advertis/DeCorr.pdf
def overCorrect(z1, z2, idx):
    defidx = np.isfinite(z1+z2)
    z11 = z1[np.logical_and(defidx, idx)]; 
    z22 = z2[np.logical_and(defidx, idx)]
    if np.sum(z1 < 0) != 0 and np.sum(z2 < 0) == 0:
        C = np.corrcoef(np.power(z11,2), z22)
        print "correlation between squred Z-score: ", C[0,1]
    elif np.sum(z1 < 0) == 0 and np.sum(z2 < 0) != 0:
        C = np.corrcoef(z11, np.power(z22,2))
        print "correlation between squred Z-score: ", C[0,1]
    else:
        C = np.corrcoef(z11, z22)
        print "correlation between Z score: ", C[0,1]
    Z = np.row_stack([z1, z2])
    print Z.shape
    z_adj = np.dot(linalg.fractional_matrix_power(C, -1/2), Z)
    print np.corrcoef(z_adj[0,:], z_adj[1,:])[0,1]
    return z2logp(z_adj[0,:]), z2logp(z_adj[1,:]) 
  
def z2logp(zscores, tails = 2):
    """
    compute coresponding -log10 p values of given z values.
    """
    return -np.log10(tails * stats.norm.cdf(-np.fabs(zscores)))
