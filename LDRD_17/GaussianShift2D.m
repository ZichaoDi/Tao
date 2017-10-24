function [y, y1,y2] = GaussianShift2D(x, delta)
% The size of the matrix.
[N, M] = size(x);
% FFT of our possibly padded input signal.
X = fft2(x);
sigma=1/2.355;

scale=1/(sqrt(2*pi)*sigma);
range_x=[0:floor(N/2)-1 floor(-N/2):-1]';
range_y=[0:floor(M/2)-1 floor(-M/2):-1];
G_x=exp(-(range_x-delta(1)).^2./(2*sigma^2));
dG_x=((range_x-delta(1))./(sigma^2)).*exp(-(range_x-delta(1)).^2./(2*sigma^2));
G_y=exp(-(range_y-delta(2)).^2./(2*sigma^2));
dG_y=((range_y-delta(2))./(sigma^2)).*exp(-(range_y-delta(2)).^2./(2*sigma^2));
y=scale^2*real(ifft2((X.*fft2(G_x*G_y))));
% Invert the FFT.
y1 = scale^2*real(ifft2((X.*fft2(dG_x*G_y))));
y2 = scale^2*real(ifft2(X.*fft2(G_x*dG_y)));

