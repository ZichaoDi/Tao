function currentSlice=align_center(currentSlice)
% load 2dSlice/Slice10;
% currentSlice=squeeze(data(40,:,:))';
% for i=1:size(currentSlice,2)-1
%     currentSlice(:,i+1)=AlignCenter(currentSlice(:,i),currentSlice(:,i+1));
% end
for i=1:size(currentSlice,1)
    currentSlice(i,:)=AlignCenter(currentSlice(i,:));
end

function alignedSignal=AlignCenter(S1)
ind=[287:390];
[m,d]=max(S1(ind));
delay=ind(d)-length(S1); 
alignedSignal=zeros(1,length(S1));
if(delay>=0)
    alignedSignal(delay+1:end)=S1(1:end-delay);
else
    alignedSignal(1:end+(delay))=S1(abs(delay)+1:end);
end

