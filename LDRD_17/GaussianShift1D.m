function [y, Daligned] = GaussianShift1D(x, shift)
global numThetan nTau 
global extNtau realInd 
% The size of the matrix.
% FFT of our possibly padded input signal.
sigma=1.6/2.355;
scale=1/(sqrt(2*pi)*sigma);
N=extNtau;
nchannel=size(x,3);
y=zeros(numThetan,nTau+1,nchannel);
Daligned=zeros(numThetan,numThetan,N,nchannel);
range_x=[0:floor(N/2)-1 floor(-N/2):-1]';
for n=1:numThetan
    G_x=exp(-(range_x-shift(n)).^2./(2*sigma^2));
    dG_x=((range_x-shift(n))./(sigma^2)).*exp(-(range_x-shift(n)).^2./(2*sigma^2));

    for ele=1:nchannel
        X = fft(squeeze(x(n,:,ele))');
        ytemp=scale*real(ifft(X.*(fft(G_x))));
        y(n,:,ele)=ytemp(realInd);
        Daligned(n,n,:,ele) = scale*real(ifft(X.*fft(dG_x)));
    end
end
Daligned=Daligned(:,:,realInd,:);
Daligned=reshape(Daligned,[numThetan,numThetan*length(realInd),nchannel]);
