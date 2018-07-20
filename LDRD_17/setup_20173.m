load 20173_data_MAPS;
slice=80;
data=squeeze(data(:,:,slice,:));
ind_i0=10;
indnan = find(isnan(squeeze(data(ind_i0,:,:))));
ind_xrt=11;
data(isnan(data))=0;
temp=squeeze(data(ind_i0,:,:));
temp(indnan)=max(temp(:));
data(ind_i0,:,:)=temp;
temp=squeeze(data(ind_xrt,:,:));
temp(indnan)=max(temp(:));
data(ind_xrt,:,:)=temp;
slice_tot = 1:9;%[1 2 3 4 9];
% slice_tot=slice_tot([2 3 5]);
NumElement=length(slice_tot);
%%%================================
data_h=[];
ang_rate=1;
tau_rate=2;
for ele=1:size(data,1)
    data_h(ele,:,:)=data(ele,1:ang_rate:end,1:tau_rate:end);
end
thetan = linspace(-30,150,size(data,2));
thetan_real = thetan(1:ang_rate:end); 
if(ndims(data_h)==2)
    data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
end
numThetan=size(data_h,2);
nTau=size(data_h,3)-1;
N = size(data_h,3);
data_xrt=data_h(ind_xrt,:,:); %% Downstream Transmission 
%%==============================================
I0=reshape(data_h(ind_i0,:,:),size(data_h,2),size(data_h,3));
% I0=repmat(max(squeeze(data_xrt),[],2),1,nTau+1);
DisR=sparse(squeeze(sum(data_xrt(:,:,:),1))');
Mt=-log(DisR./I0');
data_xrf_decom=[];
for ele=1:NumElement
    data_xrf_decom(ele,:,:)=data_h(slice_tot(ele),:,:);
end
% save('tomopy_test.mat','data_xrf_decom','data_xrt','thetan_real');

XRF_decom=permute(data_xrf_decom(:,:,:),[2 3 1]);
%%================================================
load Filter_data_spectra100;
data_xrf_raw=permute(spectra(1:tau_rate:end,1:ang_rate:end,:),[3 2 1]);
data_xrf_raw=sparse(reshape(double(data_xrf_raw),[size(data_xrf_raw,1),size(data_xrf_raw,2)*size(data_xrf_raw,3)]));
XRF_raw=data_xrf_raw';
clear spectra spectra_30_aligned data_sa data_ds data_h
iR=ones(N,N,NumElement);%load tomopy_20173
m_h=size(iR,1);
[x_ir,y_ir]=meshgrid(1:m_h);
[x_num,y_num]=meshgrid(linspace(1,m_h,N(1)));
iR_num=zeros(N(1),N(1),size(iR,3));
for ele=1:size(iR,3)
    iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,ele),x_num,y_num);
end
