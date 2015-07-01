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
mat_contents=sio.loadmat('/homes/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/ApsDataExtract/DogaSeeds/DownSampleSeeds28_elements.mat')

XTM=mat_contents['data_H']
factor=2
thetan = tomopy.angles(XTM.shape[2]/factor,0,360/factor)
XTM=XTM[:,:,:len(thetan)/factor]
XTM=np.transpose(XTM,(2,1,0))
rec_XTM = tomopy.recon(XTM,thetan, algorithm='gridrec')
# time.sleep(5.5)    # pause 5.5 seconds
print("something")
import matplotlib.pyplot as plt
# pylab.imshow(np.reshape(rec_XTM,(rec_XTM.shape[1],rec_XTM.shape[1])),cmap='gray')
slice = np.array([7, 8, 13, 16])

for i in range(0,4):
    plt.subplot(2,2,i+1)
    plt.imshow(rec_XTM[slice[i]])

plt.show()

matfile='/homes/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/ApsDataExtract/DogaSeeds/tomoRecon28_half.mat'
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
"""
