
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Convolution2D, MaxPooling2D
from keras.layers import Activation, Dropout, Flatten, Dense
import os
import numpy as np
import scipy.io as sio
from PIL import Image
from utils import nor_data, extract_patches, reconstruct_patches

import sys

topdir = r'./data/20170308_data2/top'
bottomdir = r'./data/20170308_data2/bottom'

im_train_x=np.array(Image.open(os.path.join(topdir,'2idd_0034_P_crop62x62_top.tif')));
im_train_y=np.array(Image.open(os.path.join(topdir,'2idd_0043_P_fine310x310_top.tif')));
im_test_x=np.array(Image.open(os.path.join(bottomdir,'2idd_0034_P_crop62x62_bottom.tif')));
im_test_y=np.array(Image.open(os.path.join(bottomdir,'2idd_0043_P_fine310x310_bottom.tif')));
patch_size=[15, 30]
patch_step=1
pxH,pyH=patch_size
pxh,pyh=np.multiply(5,patch_size)
nxH, nyH = im_train_x.shape
nxh, nyh = im_train_y.shape
train_x = extract_patches(im_train_x, patch_size, patch_step)
train_y = extract_patches(im_train_y, [pxh,pyh], patch_step*5)
X_train=np.reshape(train_x,(len(train_x),pxH,pyH))
Y_train=np.reshape(train_y,(len(train_y),pxh,pyh))

X_test= extract_patches(im_test_x, patch_size, patch_step)
X_test=np.reshape(X_test,(len(X_test),pxH,pyH))
Y_test= extract_patches(im_test_y, [pxh,pyh], patch_step*5)
Y_test=np.reshape(Y_test,(len(Y_test),pxh,pyh))

model = []
model = Sequential()
model.add(Dense(pxh*pyh, input_dim=pxH*pyH))
model.add(Dense(pxh*pyh, input_dim=pxH*pyH))
model.add(Dense(pxh*pyh, input_dim=pxH*pyH))
model.compile(loss='categorical_crossentropy', optimizer='sgd', metrics=['accuracy'])
b_size=1
model.fit(X_train, Y_train, nb_epoch=1, batch_size=b_size)
loss_and_metrics = model.evaluate(X_test, Y_test, batch_size=b_size)
Y_predict = model.predict(X_test,batch_size=b_size)  
Y_predict = np.reshape(Y_predict,(len(Y_predict),pxh,pyh));
img_rec = reconstruct_patches(Y_predict, [nxh,nyh], 5*patch_step)
import matplotlib.pyplot as plt   
plt.imshow(np.reshape(img_rec,(nxh,nyh)),interpolation='nearest') 
plt.show()
