import numpy as np
from utils import nor_data, extract_patches, reconstruct_patches
from models import mirror_model, test_model



def train(img_x, img_y, patch_size, patch_step, dim_img, nb_filters, nb_conv, batch_size, nb_epoch, x_test, y_test):
    """
    Function description.

    Parameters
    ----------
    parameter_01 : type
        Description.

    parameter_02 : type
        Description.

    parameter_03 : type
        Description.

    Returns
    -------
    return_01
        Description.
    """

    # img_x = nor_data(img_x)
    img_y = nor_data(img_y)
    img_input = extract_patches(img_x, patch_size, 1)
    img_output = extract_patches(img_y, patch_size, 1)
    img_input = np.reshape(img_input, (len(img_input), 1, dim_img, dim_img))
    img_output = np.reshape(img_output, (len(img_input), 1, dim_img, dim_img))

    # test_x = nor_data(x_test)
    test_y = nor_data(y_test)
    test_x = extract_patches(x_test, patch_size, 1)
    test_y = extract_patches(test_y, patch_size, 1)
    test_x = np.reshape(test_x, (len(test_x), 1, dim_img, dim_img))
    test_y = np.reshape(test_y, (len(test_y), 1, dim_img, dim_img))

    mdl = mirror_model(dim_img, nb_filters, nb_conv)
    print(mdl.summary())
    mdl.fit(img_input, img_output, batch_size=batch_size, nb_epoch=nb_epoch, validation_data=(test_x,test_y))
    return mdl


def predict(img, patch_size, patch_step, nb_filters, nb_conv, batch_size, dim_img, wpath):
    """
    the cnn model for image transformation


    Parameters
    ----------
    img : array
        The image need to be calculated

    patch_size : (int, int)
        The patches dimension

    dim_img : int
        The input image dimension

    Returns
    -------
    img_rec
        Description.

      """
    # img = nor_data(img)
    img_y, img_x = img.shape
    x_img = extract_patches(img, patch_size, patch_step)
    x_img = np.reshape(x_img, (len(x_img), 1, dim_img, dim_img))
    mdl = mirror_model(dim_img, nb_filters, nb_conv)
    mdl.load_weights(wpath)
    y_img = mdl.predict(x_img, batch_size=batch_size)
    del x_img
    y_img = np.reshape(y_img, (len(y_img), dim_img, dim_img))
    img_rec = reconstruct_patches(y_img, (img_y, img_x), patch_step)
    return img_rec
