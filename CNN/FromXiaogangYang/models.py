from keras.models import Model
from keras import backend as K
from keras.layers.core import Dense, Reshape, Flatten, Dropout, Lambda
from keras.layers import Input, merge, Conv2D, MaxPooling2D, UpSampling2D, Conv2DTranspose



def psnr(y_true, y_pred):
    return 1/(10.0 * K.log(1.0 / (K.mean(K.square(y_pred - y_true)))) / K.log(10.0))

def model_test(iy, ix, nb_filters, nb_conv):

    inputs = Input((iy, ix, 1))

    conv1 = Conv2D(nb_filters, ((nb_conv, nb_conv)), activation='relu', padding='same')(inputs)
    conv1a = Conv2D(nb_filters, (nb_conv, nb_conv),
                   strides=(2,2),activation='relu', padding='same')(conv1)



    conv2 = Conv2D(nb_filters * 2, (nb_conv, nb_conv), activation='relu', padding='same')(conv1a)
    conv2a = Conv2D(nb_filters * 2, (nb_conv, nb_conv),
                   strides=(2, 2), activation='relu', padding='same')(conv2)


    conv3 = Conv2D(nb_filters * 2, (nb_conv, nb_conv), activation='relu', padding='same')(conv2a)
    conv3a = Conv2D(nb_filters * 2, (nb_conv, nb_conv),
                   strides=(2, 2), activation='relu', padding='same')(conv3)


    conv4 = Conv2D(nb_filters * 4, (nb_conv, nb_conv), activation='relu', padding='same')(conv3a)
    conv4 = Conv2D(nb_filters * 4, (nb_conv, nb_conv), activation='relu', padding='same')(conv4)
    conv4 = Conv2D(1, (nb_conv, nb_conv), activation='relu', padding='same')(conv4)
    #
    fc1 = Flatten()(conv4)
    fc1 = Dense(ix * iy / 128, activation='relu')(fc1)
    fc1 = Dropout(0.2)(fc1)
    fc1 = Dense(ix * iy / 128, activation='relu')(fc1)
    fc1 = Dropout(0.25)(fc1)
    fc1 = Dense(ix * iy / 64, activation='relu')(fc1)
    fc1 = Dropout(0.25)(fc1)
    fc1 = Reshape((iy / 8, ix / 8, 1))(fc1)

    fc2 = Conv2DTranspose(nb_filters * 4, (nb_conv, nb_conv), activation='relu', padding='same')(fc1)

    fc2 = Conv2DTranspose(nb_filters * 8, (nb_conv, nb_conv),
                          strides=(2, 2), activation='relu', padding='same')(fc2)

    up1 = merge([fc2, conv3], mode='concat', concat_axis=3)


    conv6 = Conv2DTranspose(nb_filters * 2, (nb_conv, nb_conv), activation='relu', padding='same')(up1)
    conv6 = Conv2DTranspose(nb_filters * 2, (nb_conv, nb_conv),
                            strides=(2, 2), activation='relu', padding='same')(conv6)

    up2 = merge([conv6, conv2], mode='concat',concat_axis=3)

    conv7 = Conv2DTranspose(nb_filters * 2, (nb_conv, nb_conv), activation='relu', padding='same')(up2)
    conv7 = Conv2DTranspose(nb_filters * 2, (nb_conv, nb_conv),
                            strides=(2, 2), activation='relu', padding='same')(conv7)

    up3 = merge([conv7, conv1], mode='concat',concat_axis=3)

    conv8 = Conv2DTranspose(nb_filters, (nb_conv, nb_conv), activation='relu', padding='same')(up3)
    conv8 = Conv2DTranspose(nb_filters, (nb_conv, nb_conv), activation='relu', padding='same')(conv8)

    conv8 = Conv2DTranspose(1, (3, 3), activation='relu', padding='same')(conv8)



    mdl = Model(inputs=inputs, outputs=conv8)

    mdl.compile(loss= psnr, optimizer='Adam', metrics=['mse'])
    return mdl

def model(dim_img, nb_filters, nb_conv):

    inputs = Input((1, dim_img, dim_img))
    conv1 = MeanPadding(padding=(1,1), data_format='channels_first')(inputs)
    conv1 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv1)


    pool1 = MaxPooling2D(pool_size=(2, 2))(conv1)
    conv2 = MeanPadding()(pool1)
    conv2 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv2)

    pool2 = MaxPooling2D(pool_size=(2, 2))(conv2)
    conv3 = MeanPadding()(pool2)
    conv3 = Convolution2D(nb_filters*2, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv3)
    conv3 = MeanPadding()(conv3)
    conv3 = Convolution2D(nb_filters*2, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv3)


    #
    fc1 = Flatten()(conv3)
    fc1 = Dense((dim_img / 4) ** 2, activation='relu')(fc1)
    fc1 = Dropout(0.50)(fc1)
    fc1 = Dense((dim_img / 4) ** 2, activation='relu')(fc1)
    fc1 = Dropout(0.25)(fc1)
    fc1 = Reshape((1, dim_img / 4, dim_img / 4))(fc1)
    fc2 = MeanPadding()(fc1)
    fc2 = Deconvolution2D(nb_filters*2, nb_conv, nb_conv,
                          output_shape=(None, nb_filters*2, dim_img/4, dim_img/4),
                          activation='relu', border_mode='valid')(fc2)
    fc2 = MeanPadding()(fc2)
    fc2 = Deconvolution2D(nb_filters * 2, nb_conv, nb_conv,
                          output_shape=(None, nb_filters * 2, dim_img / 4, dim_img / 4),
                          activation='relu', border_mode='valid')(fc2)



    up1 = merge([UpSampling2D(size=(2, 2))(fc2), conv2], mode='concat', concat_axis=1)
    up1 = MeanPadding()(up1)
    conv4 = Deconvolution2D(nb_filters, nb_conv, nb_conv,
                          output_shape=(None, nb_filters, dim_img/2, dim_img/2),
                          activation='relu', border_mode='valid')(up1)


    up2 = merge([UpSampling2D(size=(2, 2))(conv4), conv1], mode='concat', concat_axis=1)
    up2 = MeanPadding()(up2)
    conv5= Deconvolution2D(nb_filters, 3, 3,
                             output_shape=(None, nb_filters, dim_img, dim_img),
                             activation='relu', border_mode='valid')(up2)

    conv5 = MeanPadding()(conv5)
    conv6 = Convolution2D(1, 3, 3, border_mode='valid')(conv5)

    mdl = Model(input=inputs, output=conv6)

    mdl.compile(loss=loss_DSSIM_theano, optimizer='Adam', metrics=[mSSIM])
    return mdl

def mirror_model(dim_img, nb_filters, nb_conv):

    inputs = Input((1, dim_img, dim_img))
    conv1 = Lambda(mirror_padding2, output_shape=padding_shape2)(inputs)
    conv1 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv1)
    conv1 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv1)

    pool1 = MaxPooling2D(pool_size=(2, 2))(conv1)
    conv2 = Lambda(mirror_padding2, output_shape=padding_shape2)(pool1)
    conv2 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv2)
    conv2 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv2)

    pool2 = MaxPooling2D(pool_size=(2, 2))(conv2)
    conv3 = Lambda(mirror_padding2, output_shape=padding_shape2)(pool2)
    conv3 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv3)
    conv3 = Convolution2D(nb_filters*2, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv3)

    #
    fc1 = Flatten()(conv3)
    fc1 = Dense((dim_img / 4) ** 2, activation='relu')(fc1)
    fc1 = Dropout(0.50)(fc1)
    fc1 = Dense((dim_img / 4) ** 2, activation='relu')(fc1)
    fc1 = Dropout(0.25)(fc1)
    fc1 = Reshape((1, dim_img / 4, dim_img / 4))(fc1)
    fc2 = Lambda(mirror_padding2, output_shape=padding_shape2)(fc1)
    fc2 = Convolution2D(nb_filters*2, nb_conv, nb_conv, activation='relu', border_mode='valid')(fc2)
    fc2 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(fc2)


    up1 = merge([UpSampling2D(size=(2, 2))(fc2), conv2], mode='concat', concat_axis=1)
    up1 = Lambda(mirror_padding2, output_shape=padding_shape2)(up1)
    conv4 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(up1)
    conv4 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv4)


    up2 = merge([UpSampling2D(size=(2, 2))(conv4), conv1], mode='concat', concat_axis=1)
    up2 = Lambda(mirror_padding2, output_shape=padding_shape2)(up2)
    conv5 = Convolution2D(nb_filters, 3, 3, activation='relu', border_mode='valid')(up2)
    conv5 = Convolution2D(nb_filters, 3, 3, activation='relu', border_mode='valid')(conv5)

    conv6 = Convolution2D(1, 1, 1, border_mode='valid')(conv5)

    mdl = Model(input=inputs, output=conv6)

    mdl.compile(loss=loss_DSSIM_theano, optimizer='Adam', metrics=[mSSIM])
    return mdl



def test_model(dim_img, nb_filters, nb_conv):

    inputs = Input((1, dim_img, dim_img))
    conv1 = Lambda(mirror_padding, output_shape=padding_shape)(inputs)
    conv1 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv1)
    conv1 = Lambda(mirror_padding, output_shape=padding_shape)(conv1)
    conv1 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv1)

    conv1 = Lambda(mirror_padding, output_shape=padding_shape)(conv1)
    conv1 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='valid')(conv1)

    # conv5 = Lambda(mirror_padding, output_shape=padding_shape)(conv1)
    conv6 = Convolution2D(1, 1, 1, border_mode='valid')(conv1)

    mdl = Model(input=inputs, output=conv6)

    mdl.compile(loss=loss_DSSIM_theano, optimizer='Adam', metrics=[mSSIM])
    return mdl


def model_l1(dim_img, nb_filters, nb_conv):

    inputs = Input((1, dim_img, dim_img))
    conv1 = Convolution2D(nb_filters, 3, 3, activation='relu', border_mode='same')(inputs)
    conv1 = Convolution2D(nb_filters, 5, 5, activation='relu', border_mode='same')(conv1)

    conv5 = Deconvolution2D(nb_filters, 5, 5,
                             output_shape=(None, nb_filters, dim_img, dim_img),
                             activation='relu', border_mode='same')(conv1)
    conv5 = Deconvolution2D(nb_filters, 3, 3,
                             output_shape=(None, nb_filters, dim_img, dim_img),
                             activation='relu', border_mode='same')(conv5)
    conv6 = Convolution2D(1, 3, 3, border_mode='same')(conv5)

    mdl = Model(input=inputs, output=conv6)

    mdl.compile(loss='mean_squared_error', optimizer='Adam')
    return mdl

def model64(dim_img, nb_filters, nb_conv):

    inputs = Input((1, dim_img, dim_img))
    conv1 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='same')(inputs)
    conv1 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='same')(conv1)

    pool1 = MaxPooling2D(pool_size=(2, 2))(conv1)
    conv2 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='same')(pool1)

    pool2 = MaxPooling2D(pool_size=(2, 2))(conv2)
    conv3 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='same')(pool2)

    pool3 = MaxPooling2D(pool_size=(2, 2))(conv3)
    conv4 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='same')(pool3)

    #
    fc1 = Flatten()(conv4)
    fc1 = Dense((dim_img / 8) ** 2, activation='relu')(fc1)
    fc1 = Dropout(0.25)(fc1)
    fc1 = Dense((dim_img / 8) ** 2, activation='relu')(fc1)
    fc1 = Dropout(0.25)(fc1)
    fc1 = Reshape((1, dim_img / 8, dim_img / 8))(fc1)
    fc2 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='same')(fc1)

    up1 = merge([UpSampling2D(size=(2, 2))(fc2), conv3], mode='concat', concat_axis=1)
    conv5 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='same')(up1)

    up2 = merge([UpSampling2D(size=(2, 2))(conv5), conv2], mode='concat', concat_axis=1)
    conv6 = Convolution2D(nb_filters, nb_conv, nb_conv, activation='relu', border_mode='same')(up2)

    up3 = UpSampling2D(size=(2, 2))(conv6)
    conv6a = Convolution2D(nb_filters/2, nb_conv, nb_conv, activation='relu', border_mode='same')(up3)
    conv6b = Convolution2D(nb_filters/2, 5, 5, activation='relu', border_mode='same')(up3)
    conv6c = Convolution2D(nb_filters/2, 7, 7, activation='relu', border_mode='same')(up3)
    mg1 = merge([conv6a, conv6b, conv6c, conv1], mode='concat', concat_axis=1)

    conv7 = Convolution2D(1, 1, 1, border_mode='same')(mg1)

    mdl = Model(input=inputs, output=conv7)

    mdl.compile(loss='mean_squared_error', optimizer='Adam')
    return mdl
