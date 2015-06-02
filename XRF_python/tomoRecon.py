# _*_ coding: utf-8 _*_
"""
Created on Mon May 11 12:30:02 2015

@author: Wendydi
"""

import tomopy
import scipy.io as sio
import numpy as np

# mat_contents=sio.loadmat('/homes/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/2xfm1211_14/3dSet1741.mat')
# thetan=mat_contents['thetan']
mat_contents=sio.loadmat('/homes/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/ApsDataExtract/DogaSeeds/DownSampleSeeds111_elements.mat')

XTM=mat_contents['data_H']
thetan = tomopy.angles(XTM.shape[2]/2,0,180)
XTM=XTM[:,7:8,:len(thetan)/2]
XTM=np.transpose(XTM,(2,1,0))
rec_XTM = tomopy.recon(XTM,thetan, algorithm='gridrec')

import pylab
# pylab.imshow(np.reshape(rec_XTM,(rec_XTM.shape[1],rec_XTM.shape[1])),cmap='gray')
pylab.imshow(rec_XTM[0])
pylab.show()
# matfile='/Users/Wendydi/Tao/XRF_XTM_Simulation/data/2xfm1211_14/tomoRecon1741.mat'
# sio.savemat(matfile, mdict={'out': rec_XTM}, oned_as='row')
"""
XRF=mat_contents['XRF']
n_elements = XRF.shape[2]
rec_XRF=np.ones((51, 1741, 1741,n_elements))
for e in range(0,n_elements):
    XF=np.transpose(np.squeeze(XRF[:,:,e,:]),axes=(2,1,0))
    rec_XRF[:,:,:,e] = tomopy.recon(XF,thetan, algorithm='art')    
matfile='/homes/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/2xfm1211_14/tomoRecon1741_XRF.mat'
sio.savemat(matfile, mdict={'out': rec_XRF}, oned_as='row')
"""
