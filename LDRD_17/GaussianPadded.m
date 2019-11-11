function output=GaussianPadded(x,padwidthx,padwidthy)
global numThetan nTau nchannel NumElement
global realIndSlice realInd extNtau extNslice
sigma=1.5/2.355;
scale=1/(sqrt(2*pi)*sigma);
nstart=1*floor(padwidthx/2);
extNslice=nchannel+padwidthy;
realIndSlice=1:nchannel;
% realIndSlice=nstart+1:nchannel+nstart;
padIndSlice=setdiff([1:extNslice],realIndSlice);
extNtau=nTau+1+padwidthx;
realInd=nstart+1:nTau+1+nstart;
padInd=setdiff([1:extNtau],realInd);
N=extNtau;
M=extNslice;
for n=1:numThetan
    temp=mean(mean(x(n,:,:)))*ones(extNtau,extNslice);
    temp(realInd,realIndSlice)=x(n,:,:);
    temp1=gaussfilt2D(temp,50*[sigma,sigma]);
    temp1(realInd,realIndSlice)=x(n,:,:);
    output(n,:,:)=temp1;
end

