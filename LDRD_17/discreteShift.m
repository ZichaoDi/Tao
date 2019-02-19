for slice=1:NumElement
    tt=Mt(:,:,slice);
    alignedDiscrete=zeros(numThetan,nTau+1);
for i = 1:numThetan
    delay=round(shiftx(i));
    %%---------------- Continuous: Gaussian Convolution
    sigma=1.5/2.355;
    nT=nTau+1;
    %%---------------- Discrete: Permutation Matrix
    if(shiftx(i)>=0)
        alignedDiscrete(i,delay+1:end)=tt(i,1:end-delay);
        alignedDiscrete(i,1:delay)=tt(i,end-delay+1:end);
    else
        alignedDiscrete(i,1:end+delay)=tt(i,-delay+1:end);
        alignedDiscrete(i,end+delay+1:end)=tt(i,1:-delay);
    end
    %%--------------------------------------------------
end
aligned(:,:,slice)=alignedDiscrete;
end
