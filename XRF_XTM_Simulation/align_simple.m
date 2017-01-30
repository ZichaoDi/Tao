% load spectra_30
% spectra_30_aligned=spectra_30;
load Slice30
data(isnan(data))=0;
data(isinf(data))=0;
data=permute(data,[3 2 1]);
% load per_shift.mat
per_shift=r;
data_aligned=data;
for ch=1:size(data,3)
    currentSlice=squeeze(data(:,:,ch));
    for i=1:size(currentSlice,2)-1
        delay=-per_shift(i); 
        alignedSignal=currentSlice(:,i+1);
        if(delay>=0)
            alignedSignal(delay+1:end)=alignedSignal(1:end-delay);
        else
            alignedSignal(1:end+(delay))=alignedSignal(abs(delay)+1:end);
        end
        currentSlice(:,i+1)=alignedSignal;
    end
    data_aligned(:,:,ch)=currentSlice;
end

return;
b=zeros(73,51,1748);
for i=1:73,b(i,:,:)=imread(['2xfm_0',num2str(i+172),'_x_coord.tif']);end 
b=squeeze(b(:,30,:));
shift=zeros(numThetan,1);
for i=1:numThetan; 
    [co,la]=max(b(i,:));
    if(la<length(b(i,:))/2),shift(i)=-la;
    else,shift(i)=1748-la;end
end
