function output=GaussianPadded(x)
global numThetan nTau nslice NumElement
global realIndSlice realInd extNtau extNslice
sigma=1.5/2.355;
scale=1/(sqrt(2*pi)*sigma);
padwidth=0;%1*max(nTau+1,nslice);
nstart=1*floor(padwidth/2);
extNslice=nslice+padwidth;
realIndSlice=nstart+1:nslice+nstart;
padIndSlice=setdiff([1:extNslice],realIndSlice);
extNtau=nTau+1+padwidth;
realInd=nstart+1:nTau+1+nstart;
padInd=setdiff([1:extNtau],realInd);
N=extNtau;
M=extNslice;
for n=1:numThetan
    for ele=1:NumElement
        temp=zeros(extNtau,extNslice);temp(realInd,realIndSlice)=x(n,:,:,ele);
        temp1=gaussfilt2D(temp,5.*[sigma,sigma]);
        temp(padInd,padIndSlice)=temp1(padInd,padIndSlice); 
        output(n,:,:,ele)=temp;
    end
end
