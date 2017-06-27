load ./data/ApsDataExtract/DogaSeeds/DownSampleSeeds441_elements.mat
load RawSpectra_seed.mat
if(angleScale==1)
    cutInd=floor(size(data_H,3)/2);
    data_H=data_H(:,:,1:cutInd);
    spectra=spectra(:,:,1:cutInd);
end
data = permute(data_H,[2,3,1]);
clear data_H
slice =[9 13 14 17];%5:27;%      
data_h=[];
ang_rate=2;
tau_rate=20;
for ele=1:size(data,1)
    data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
end
if(ndims(data_h)==2)
    data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
end
data_xrt=squeeze(data_h(3,:,:)); %transimission data
% I0=max(data_xrt(:));
I0=reshape(data_h(2,:,:),size(data_h,2),size(data_h,3));% 
% DisR=sparse(data_xrt');
load DisR_removeStr
data_xrf_decom=data_h(slice,:,:);
numThetan=size(data_h,2);
nTau=size(data_h,3)-1;
N = size(data_h,3);
data_xrf_raw=permute(spectra(1:tau_rate:end,:,1:ang_rate:end),[2 3 1]);
data_xrf_raw=sparse(reshape(double(data_xrf_raw),[size(data_xrf_raw,1),size(data_xrf_raw,2)*size(data_xrf_raw,3)]));
[x_ir,y_ir]=meshgrid(1:size(iR,1));
[x_num,y_num]=meshgrid(linspace(1,size(iR,1),N(1)));
iR_num=zeros(N(1),N(1),length(slice));
for ele=1:size(iR_num,3)
    iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,slice(ele)),x_num,y_num);
end
save('tomopytest.mat','data_xrf_decom');
