load Paunesku_tomo86;
ind_i0=2;
ind_xrt=3;
slice_tot = [7 8 11 14];
%%%================================
load tomoRod
data_h=[];
ang_rate=1;
tau_rate=1;
thetan=ang(timeInd);
data=data(:,timeInd,:);
for ele=1:size(data,1)
    data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
    % data_h(ele,:,:)=sum_interval(squeeze(data(ele,:,:)),ang_rate,tau_rate);
end
thetan_real = thetan(1:ang_rate:end)'; 
if(ndims(data_h)==2)
    data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
end
numThetan=size(data_h,2);
nTau=size(data_h,3)-1;
N = size(data_h,3);
data_xrt=data_h(ind_xrt,:,:); %% Downstream Transmission 
%%==============================================
I0=repmat(max(squeeze(data_xrt),[],2),1,nTau+1);
DisR=sparse(squeeze(sum(data_xrt(:,:,:),1))');
data_xrf_decom=[];
for ele=1:length(slice_tot)
data_xrf_decom(ele,:,:)=data_h(slice_tot(ele),:,:);
end
% save('tomopy_test.mat','data_xrf_decom','data_xrt','ang');
% return;
XRF_decom=permute(data_xrf_decom(:,:,:),[2 3 1]);
%%================================================
load tomopy_paunesku;
m_h=size(iR,1);
[x_ir,y_ir]=meshgrid(1:m_h);
[x_num,y_num]=meshgrid(linspace(1,m_h,N(1)));
iR_num=zeros(N(1),N(1),size(iR,3));
for ele=1:size(iR,3)
    iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,ele),x_num,y_num);
end
