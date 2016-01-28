
# coding: utf-8

# In[11]:
from PIL import Image
import os
import tomopy
import numpy as np
import scipy.io as sio
"""
fname = '/nfs2/wendydi/Documents/Research/APSdata/Seed/sample_name.h5'
data = tomopy.read_hdf5(fname, '/xrfmap/roimap/sum_cor')
"""
data=[]
for fn in os.listdir('./'):
    im = Image.open(fn)
    data.append(np.array(im))
data = np.array(data)
print data.shape, data.max(), data.min()

import pylab
# get_ipython().magic(u'matplotlib inline')
# data = data[0:360,:,:]
# data = data[0::25,:,:]
# data = np.transpose(data,(0,2,1))
num_projections, num_elements, num_pixels = data.shape
# num_projections = 24
ang = tomopy.angles(num_projections, 0, 180)
center=np.array([100,100])
rec = tomopy.recon(data[0:len(ang), :, :], ang,center,emission=True, algorithm='gridrec')

import matplotlib.pyplot as plt
slice = np.array([7, 8, 13, 16])
for i in range(0,6):
    plt.subplot(2,3,i+1)
    plt.imshow(rec[i])

plt.show()

"""
fname = '/home/oxygen/DGURSOY/Data/Trunk/aps13id/xfm_data_raw.h5'
raw_data = tomopy.read_hdf5(fname, 'exchange/data', slc=((10, 11), None))
print raw_data.shape, raw_data.max(), raw_data.min()


# In[31]:

import pylab
# get_ipython().magic(u'matplotlib inline')

pylab.plot(np.log(raw_data[0].sum(axis=1)))
pylab.grid('on')
pylab.title('Spectrum of a region')
pylab.xlabel('Energy')
pylab.ylabel('Total counts')
pylab.show()


# In[ ]:
"""


