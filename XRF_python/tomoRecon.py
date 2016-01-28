# _*_ coding: utf-8 _*_
"""
Created on Mon May 11 12:30:02 2015

@author: Wendydi
"""

import tomopy
import os
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
# mat_contents=sio.loadmat('/homes/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/2xfm1211_14/3dSet1741.mat')
# thetan=mat_contents['thetan']
mat_contents=sio.loadmat('/homes/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/ApsDataExtract/DogaSeeds/DownSampleSeeds221_elements.mat')

center=np.array([842])
rec_XTM=[]
for fn in range(0,1): # os.listdir('./'):
    # mat_contents=sio.loadmat('Slice1_59')
    XTM = np.array(mat_contents['data_H'])
    XTM.shape
    factor=1
    # XTM=np.expand_dims(XTM,axis=2)
    XTM=np.transpose(XTM,(1,0,2))
    thetan = tomopy.angles(XTM.shape[0],0,360/factor)
    # rec_XTM.append(tomopy.recon(XTM,thetan,algorithm='mlem',num_iter=10))
    rec_XTM=tomopy.recon(XTM,thetan,algorithm='mlem',num_iter=10)

rec_XTM=np.array(rec_XTM)
rec_XTM.shape
# pylab.imshow(np.reshape(rec_XTM,(rec_XTM.shape[1],rec_XTM.shape[1])),cmap='gray')
slice = np.arange(0,6)
for i in range(0,45):
    plt.subplot(5,9,i+1)
    plt.imshow(np.squeeze(rec_XTM[i,:,:]))

plt.show()

"""
matfile='/homes/wendydi/Documents/Research/APSdata/GlassRod/2dSlice/rec_XTM_59.mat'
sio.savemat(matfile, mdict={'out': rec_XTM}, oned_as='row')
XRF=mat_contents['XRF']
n_elements = XRF.shape[2]
rec_XRF=np.ones((51, 1741, 1741,n_elements))
for e in range(0,n_elements):
    XF=np.transpose(np.squeeze(XRF[:,:,e,:]),axes=(2,1,0))
    rec_XRF[:,:,:,e] = tomopy.recon(XF,thetan, algorithm='art')    
matfile='/homes/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/2xfm1211_14/tomoRecon1741_XRF.mat'
sio.savemat(matfile, mdict={'out': rec_XRF}, oned_as='row')
"""
