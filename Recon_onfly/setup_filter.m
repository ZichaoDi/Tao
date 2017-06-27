load Filter_data100;
ind_i0=29;
indnan = find(isnan(squeeze(data(ind_i0,:,:))));
ind_xrt=30;
data(isnan(data))=0;
        temp=squeeze(data(ind_i0,:,:));
        temp(indnan)=max(temp(:));
        data(ind_i0,:,:)=temp;
        temp=squeeze(data(ind_xrt,:,:));
        temp(indnan)=max(temp(:));
        data(ind_xrt,:,:)=temp;
slice_tot = [2 3 4 8 16 19];
%%%================================
data_h=[];
ang_rate=1;
tau_rate=6;
for ele=1:size(data,1)
    data_h(ele,:,:)=data(ele,1:ang_rate:end,1:tau_rate:end);
end
thetan = linspace(1,360,size(data,2));
thetan_real = thetan(1:ang_rate:end); 
if(ndims(data_h)==2)
    data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
end
numThetan=size(data_h,2);
thetan=linspace(0,360,numThetan);
nTau=size(data_h,3)-1;
N = size(data_h,3);
data_xrt=data_h(ind_xrt,:,:); %% Downstream Transmission 
%%==============================================
I0=reshape(data_h(ind_i0,:,:),size(data_h,2),size(data_h,3));
DisR=sparse(squeeze(sum(data_xrt(:,:,:),1))');
data_xrf_decom=[];
for ele=1:NumElement
    data_xrf_decom(ele,:,:)=data_h(slice_tot(ele),:,:);
end
% save('tomopy_test.mat','data_xrf_decom','data_xrt');
XRF_decom=permute(data_xrf_decom(:,:,:),[2 3 1]);
%%================================================
load Filter_data_spectra100;
data_xrf_raw=permute(spectra(1:tau_rate:end,1:ang_rate:end,:),[3 2 1]);
data_xrf_raw=sparse(reshape(double(data_xrf_raw),[size(data_xrf_raw,1),size(data_xrf_raw,2)*size(data_xrf_raw,3)]));
XRF_raw=data_xrf_raw';
clear spectra spectra_30_aligned data_sa data_ds data_h
load tomopy_filter
m_h=size(iR,1);
[x_ir,y_ir]=meshgrid(1:m_h);
[x_num,y_num]=meshgrid(linspace(1,m_h,N(1)));
iR_num=zeros(N(1),N(1),size(iR,3));
for ele=1:size(iR,3)
    iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,ele),x_num,y_num);
end
