
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Convolution2D, MaxPooling2D
from keras.layers import Activation, Dropout, Flatten, Dense
import os
import numpy as np
import scipy.io as sio

train=sio.loadmat('ML_train_36000.mat')
X_train=train['sino']
Y_train=train['W']

X_train = np.transpose(X_train)
Y_train = np.transpose(Y_train)

test=sio.loadmat('ML_test.mat')
X_test=test['sino']
Y_test=test['W']

X_test = np.transpose(X_test)
Y_test = np.transpose(Y_test)

model = []
model = Sequential()
model.add(Dense(2500, input_dim=1500))
model.compile(loss='categorical_crossentropy', optimizer='sgd', metrics=['accuracy'])
model.fit(X_train, Y_train, nb_epoch=100, batch_size=32)
loss_and_metrics = model.evaluate(X_test, Y_test, batch_size=32)
Y_predict = model.predict(np.reshape(X_train[7846],(1,1500)),batch_size=1)  
import matplotlib.pyplot as plt   
plt.imshow(np.reshape(Y_predict,(50,50)),interpolation='nearest') 
plt.show()
