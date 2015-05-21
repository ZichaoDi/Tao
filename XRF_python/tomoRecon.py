# _*_ coding: utf-8 _*_
"""
Created on Mon May 11 12:30:02 2015

@author: Wendydi
"""

import tomopy
import scipy.io as sio
import numpy as np

mat_contents=sio.loadmat('/homes/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/2xfm1211_14/3dSet1741.mat')
thetan=mat_contents['thetan']
"""
XTM=mat_contents['XTM']

XTM=np.transpose(XTM,axes=(2,1,0))
rec_XTM = tomopy.recon(XTM,thetan, algorithm='art')

# import pylab
# pylab.imshow(rec[64],cmap='gray')

matfile='/Users/Wendydi/Tao/XRF_XTM_Simulation/data/2xfm1211_14/tomoRecon1741.mat'
sio.savemat(matfile, mdict={'out': rec_XTM}, oned_as='row')
"""
XRF=mat_contents['XRF']
n_elements = XRF.shape[2]
rec_XRF=np.ones((51, 1741, 1741,n_elements))
for e in range(0,n_elements):
    XF=np.transpose(np.squeeze(XRF[:,:,e,:]),axes=(2,1,0))
    rec_XRF[:,:,:,e] = tomopy.recon(XF,thetan, algorithm='art')
    
matfile='/homes/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/2xfm1211_14/tomoRecon1741_XRF.mat'
sio.savemat(matfile, mdict={'out': rec_XRF}, oned_as='row')
