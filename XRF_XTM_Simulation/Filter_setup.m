load data/ApsDataExtract/Luxi/Filter_data100;
ind_i0=43;
ind_xrt=38;
slice_tot = [3 4 25 30]; %GlassRod%
% load slice30
% slice_tot = [30 31 20 29]; %GlassRod%
% ind_i0=2;
% ind_xrt=3;
% data(isinf(data))=0;
% data(isnan(data))=0;
% xrt=-log(data(ind_xrt,:,:)./data(ind_i0,:,:));
% xrf(1,:,:)=sum(data([slice_tot(1) slice_tot(2)],:,:),1);
% xrf(2,:,:)=data(slice_tot(3),:,:);
% xrf(3,:,:)=data(slice_tot(4),:,:);
% save('tomopytest.mat','xrf','xrt');
load tomoRod
data_h=[];
ang_rate=2;
tau_rate=25;
for ele=1:size(data,1)
    data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
end
thetan = linspace(1,360,size(data,2));
thetan_real = thetan(1:ang_rate:end); 
if(ndims(data_h)==2)
    data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
end
numThetan=size(data_h,2);
nTau=size(data_h,3)-1;
N = size(data_h,3);
I0=reshape(data_h(ind_i0,:,:),size(data_h,2),size(data_h,3));
data_xrt=data_h(ind_xrt,:,:); %% Downstream Transmission 
DisR=sparse(squeeze(sum(data_xrt(:,:,:),1))');
data_xrf_decom=[];
data_xrf_decom(1,:,:)=data_h(slice_tot(1),:,:)+data_h(slice_tot(2),:,:);
data_xrf_decom(2,:,:)=data_h(slice_tot(3),:,:);
data_xrf_decom(3,:,:)=data_h(slice_tot(4),:,:);
% save('tomopy_coarse.mat','data_xrf_decom','data_xrt');
% return;
XRF_decom=permute(data_xrf_decom(:,:,:),[2 3 1]);
load spectra_30_aligned;
data_xrf_raw=permute(spectra_30_aligned(1:tau_rate:end,1:ang_rate:end,:),[3 2 1]);
data_xrf_raw=sparse(reshape(double(data_xrf_raw),[size(data_xrf_raw,1),size(data_xrf_raw,2)*size(data_xrf_raw,3)]));
XRF_raw=data_xrf_raw';
clear spectra_30_aligned data_sa data_ds data_h
m_h=size(iR,1);
[x_ir,y_ir]=meshgrid(1:m_h);
[x_num,y_num]=meshgrid(linspace(1,m_h,N(1)));
iR_num=zeros(N(1),N(1),size(iR,3));
for ele=1:size(iR,3)
    iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,ele),x_num,y_num);
end

