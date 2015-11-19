function currentSlice=align_center(currentSlice)
% load 2dSlice/Slice10;
% currentSlice=squeeze(data(40,:,:))';
for i=1:size(currentSlice,2)-1
    currentSlice(:,i+1)=AlignCenter(currentSlice(:,i),currentSlice(:,i+1));
end

    function alignedSignal=AlignCenter(S1,S2)
    [cor,lag]=xcorr(S1,S2);
    [m,d]=max(cor); 
    delay=d-max(length(S1),length(S2)); 
    alignedSignal=zeros(size(S2));
    if(delay>=0)
        alignedSignal(delay+1:end)=S2(1:end-delay);
    else
        alignedSignal(1:end+(delay))=S2(abs(delay)+1:end);
    end

