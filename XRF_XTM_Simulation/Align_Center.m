function shift1=Align_Center(currentSlice)
% load Slice10;
% currentSlice=squeeze(data(40,:,:))';
for i=1:size(currentSlice,2)-1
    %currentSlice(:,i+1)=AlignCenter(currentSlice(:,i),currentSlice(:,i+1));
    shift1(i)=AlignCenter(currentSlice(:,1),currentSlice(:,i+1));
end





function delay=AlignCenter(S1,S2)

[cor,lag]=xcorr(S1,S2);
[m,d]=max(cor); 
delay=d-max(length(S1),length(S2)); 
alignedSignal=zeros(1,length(S1));
if(delay>=0)
    alignedSignal(delay+1:end)=S1(1:end-delay);
else
    alignedSignal(1:end+(delay))=S1(abs(delay)+1:end);
end
