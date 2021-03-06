function [y, Daligned] = GaussianShift2D(x, shift)
global numThetan nTau nslice NumElement
global extNtau extNslice realInd realIndSlice
y=zeros(numThetan,nTau+1,nslice,NumElement);
% Daligned=zeros(2*numThetan,numThetan*(nTau+1),nslice,NumElement);
% The size of the matrix.
delta=reshape(shift,2,numThetan);
% FFT of our possibly padded input signal.
sigma=1.5/2.355;
scale=1/(sqrt(2*pi)*sigma);
N=extNtau;
M=extNslice;
Daligned=zeros(2*numThetan,numThetan,N,M,NumElement);
range_x=[0:floor(N/2)-1 floor(-N/2):-1]';
range_y=[0:floor(M/2)-1 floor(-M/2):-1];
for n=1:numThetan
    G_x=exp(-(range_x-delta(1,n)).^2./(2*sigma^2));
    dG_x=((range_x-delta(1,n))./(sigma^2)).*exp(-(range_x-delta(1,n)).^2./(2*sigma^2));
    G_y=exp(-(range_y-delta(2,n)).^2./(2*sigma^2));
    dG_y=((range_y-delta(2,n))./(sigma^2)).*exp(-(range_y-delta(2,n)).^2./(2*sigma^2));

    for ele=1:NumElement
        X = fft2(squeeze(x(n,:,:,ele)));
        ytemp=scale^2*real(ifft2(X.*(fft2(G_x*G_y))));
        y(n,:,:,ele)=ytemp(realInd,realIndSlice);
        Daligned(2*n-1,n,:,:,ele) = scale^2*real(ifft2(X.*fft2(dG_x*G_y)));
        Daligned(2*n,n,:,:,ele) = scale^2*real(ifft2(X.*fft2(G_x*dG_y)));
    end
end
Daligned=Daligned(:,:,realInd,realIndSlice,:);
