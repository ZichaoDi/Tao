function currentSlice=align_center(currentSlice,shift1)
global i
% load 2dSlice/Slice10;
% currentSlice=squeeze(data(40,:,:))';
for i=1:size(currentSlice,2)-1
    currentSlice(:,i+1)=AlignCenter(currentSlice(:,i),currentSlice(:,i+1),shift1);
end

function alignedSignal=AlignCenter(S1,S2,shift1)
global i
[cor,lag]=xcorr(S1,S2);
[m,d]=max(cor);
% [peaks,loc]=findpeaks(double(cor));
% sub_ind=find(loc>1720 & loc<1750);
% [m,d]=max(peaks(sub_ind)); 
%  figure(3);plot(cor);
% hold on; plot(loc(sub_ind),peaks(sub_ind),'ro')
% pause;
% % hold off; axis([0 1750 0 6e4])
% d=loc(sub_ind(d));
delay=-shift1(i);% d-max(length(S1),length(S2)); 
alignedSignal=S2;
if(delay>=0)
    alignedSignal(delay+1:end)=S2(1:end-delay);
else
    alignedSignal(1:end+(delay))=S2(abs(delay)+1:end);
end

