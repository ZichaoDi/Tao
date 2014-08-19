%============= Visulize hdf and tif data
close all;
warning off;
T1=double(imread('zeolite_4A_2_5x_00091.tif'));
T2=double(imread('zeolite_4A_2_5x_00092.tif'));
Tc1=GrayScale(T1);
Tc2=GrayScale(T2);
T=T2./T1;
Tc=GrayScale(T);
figure, imshow(uint8(Tc1)); figure, imshow(uint8(Tc2));
figure,
imshow(uint8(Tc))