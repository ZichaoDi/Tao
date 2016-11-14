function [nima]=xRemoveStripes(ima,sigma)

 % wavelet decomposition
[af,sf]=farras;
w=dwt2D(ima,1,af);

% FFT transform of horizontal frequency bands
for i=1:3
fC=fftshift(fft(w{1}{i}));
[my,mx]=size(fC);

% damping of vertical stripe information
damp=1-exp(-[-floor(my/2):-floor(my/2)+my-1].^2/( 2 * sigma^2));
fC=fC.* repmat(damp',1,mx);
% inverse FFT
w{1}{i}=ifft(ifftshift(fC));
end

% wavelet reconstruction
nima=idwt2D(w,1,sf);
