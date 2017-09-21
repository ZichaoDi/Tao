import dxchange
import numpy as np
from transform import train, predict


batch_size = 400
nb_epoch = 50
dim_img = 64
nb_filters = 40
nb_conv = 3
patch_step = 1
patch_size = (dim_img, dim_img)
smooth = 1.

test1 = dxchange.read_tiff('/home/beams/YANGX/ptychography/datatolearn/lr_0034.tiff')
ih, iw = test1.shape
img_prd = np.zeros((ih, iw))
wpath = 'th_weights/full2.h5'
img1 = test1[:ih/2, :iw/2]
img_rec1 = predict(img1, patch_size, patch_step, nb_filters, nb_conv, batch_size, dim_img, wpath)
img_prd[:ih/2, :iw/2] = img_rec1

img1 = test1[ih/2+1:, iw/2+1:]
img_rec1 = predict(img1, patch_size, patch_step, nb_filters, nb_conv, batch_size, dim_img, wpath)
img_prd[ih/2+1:, iw/2+1:] = img_rec1

img1 = test1[:ih/2, iw/2+1:]
img_rec1 = predict(img1, patch_size, patch_step, nb_filters, nb_conv, batch_size, dim_img, wpath)
img_prd[:ih/2, iw/2+1:] = img_rec1

img1 = test1[ih/2+1:, :iw/2]
img_rec1 = predict(img1, patch_size, patch_step, nb_filters, nb_conv, batch_size, dim_img, wpath)
img_prd[ih/2+1:, :iw/2] = img_rec1
fsave = '/home/beams/YANGX/ptychography/full0034_prd'
dxchange.write_tiff(img_prd, fsave, dtype='float32')
