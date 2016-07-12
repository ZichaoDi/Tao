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
% alignedSignal=zeros(size(S2));
% alignedSignal(delay+1:end)=S2(1:end-delay);
