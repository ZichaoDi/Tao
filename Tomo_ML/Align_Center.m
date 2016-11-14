function [shift,aligned]=Align_Center(currentSlice)
aligned(:,1)=currentSlice(:,1);
for i=1:size(currentSlice,2)-1
    [shift(i),aligned(:,i+1)]=AlignCenter(currentSlice(:,i),currentSlice(:,i+1));
end

function [delay,alignedSignal]=AlignCenter(S1,S2)

[cor,lag]=xcorr(S1,S2);
[~,d]=max(cor); 
delay=lag(d);%d-max(length(S1),length(S2)); 
alignedSignal=zeros(1,length(S1));
if(delay>=0)
    alignedSignal(delay+1:end)=S1(1:end-delay);
else
    alignedSignal(1:end+(delay))=S1(abs(delay)+1:end);
end
