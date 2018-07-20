alignedDiscrete=zeros(nTau+1,numThetan,nslice);
for i = 1:numThetan
    delay=round(shiftx(i));
    %%---------------- Continuous: Gaussian Convolution
    sigma=1.5/2.355;
    nT=nTau+1;
    %%---------------- Discrete: Permutation Matrix
    if(delay>=0)
        alignedDiscrete(delay+1:end,i,slice)=squeeze(Mt(i,slice,1:end-delay))';
        alignedDiscrete(slice,1:delay,i,slice)=squeeze(Mt(i,slice,end-delay+1:end))';
    else
        alignedDiscrete(1:end+delay,i,slice)=squeeze(Mt(i,slice,-delay+1:end))';
        alignedDiscrete(end+delay+1:end,i,slice)=squeeze(Mt(i,slice,1:-delay))';
    end
    %%--------------------------------------------------
end
