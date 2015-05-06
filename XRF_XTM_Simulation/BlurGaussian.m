function blurimage=BlurGaussian(origimage)
%Blur Kernel

[h, w] = size(origimage);

%Gaussian Blur
s = 1; %variance
m = [w/2,h/2]; % mean
[X, Y] = meshgrid(1:w,1:h);
kernel = (1/(2*pi*s^2))*exp(-((X-m(1)).^2 + (Y-m(2)).^2)/(2*s^2));

kernelimage = kernel;%zeros(h,w);
% kernelimage(1:ksize, 1:ksize) = kernel;

%Perform 2D FFTs
fftimage = fft2(double(origimage));
fftkernel = fft2(kernelimage);

%Set all zero values to minimum value
fftkernel(fftkernel == 0) = 1e-6;

%Multiply FFTs
fftblurimage = fftimage.*fftkernel;

%Perform Inverse 2D FFT
blurimage=rot90(ifft2(fftblurimage),2);

close all;
figure, imagesc(origimage)
axis square
colormap gray
title('Original Image')
set(gca, 'XTick', [], 'YTick', [])
%Display Kernel
figure, imagesc(kernel)
axis square
title('Blur Kernel')
colormap gray
%Display Blurred Image
figure, imagesc(blurimage)
axis square
title('Blurred Image')
colormap gray
set(gca, 'XTick', [], 'YTick', [])