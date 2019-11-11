function data_aligned=align_simple(data,shift)
%%=============data: nTau x numThetan x nchannel
% load slice30
% data(isnan(data))=0;
% data(isinf(data))=0;
% data=permute(data,[3 2 1]);
% load per_shift.mat
data_aligned=data;
nTau=size(data,1)-1;
sigma=1.5/2.355;
scale=1/(sqrt(2*pi)*sigma);
range=[0:floor((nTau+1)/2)-1 floor(-(nTau+1)/2):-1]';% / (nTau+1);
discrete=0;
for ch=1:size(data,3)
    currentSlice=squeeze(data(:,:,ch));
    for i=1:size(currentSlice,2)
        delay=shift(i); 
        alignedSignal=currentSlice(:,i);
        if(discrete)
            if(delay>=0)
                alignedSignal(delay+1:end)=alignedSignal(1:end-delay);
            else
                alignedSignal(1:end+(delay))=alignedSignal(abs(delay)+1:end);
            end
        else
            G=exp(-(range-shift(i)).^2./(2*sigma^2));
            alignedSignal=scale*real(ifft((fft(G)).*(fft(alignedSignal))));
        end
        currentSlice(:,i)=alignedSignal;
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
