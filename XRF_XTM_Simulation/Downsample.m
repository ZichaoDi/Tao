% function [iR]=Downsample(rate,spectra)
load RawSpectra
rate=1;
%% Downsample the experimental data
%% sample rate: 2n-1
data=h5read('../../APSdata/Seed/xfm_data_elements.h5','/exchange/data');
data_H=[];
%rate=15;
ang=linspace(0,360,725);  
%for i=1:size(data,2) % Number of elements
%for j=1:size(data,3) % Number of angles
% data_H(:,i,j)=downsample(data(:,i,j),rate);
%end
%end

spec_H=[];
for i=1:size(spectra,1)
if(i*rate<=size(spectra,1))
spec_H(i,:,:)=sum(spectra(rate*(i-1)+1:i*rate,:,:),1);
data_H(i,:,:)=sum(data(rate*(i-1)+1:i*rate,:,:),1);
else
spec_H(i,:,:)=sum(spectra(rate*(i-1)+1:end,:,:),1);
data_H(i,:,:)=sum(data(rate*(i-1)+1:end,:,:),1);
break;
end
end

iR=[];for i=1:27,iR(:,:,i)=iradon(squeeze(data_H(:,i,:)),ang);end
save(['./data/ApsDataExtract/DogaSeeds/DownSampleSeeds',num2str(size(data_H,1)),'_elements.mat'],'data_H','iR','spec_H'); 
