load(['miller_raw_fourier',num2str(slice),'.mat']);
load ang_miller;
nslice=74;
off_ind=setdiff(1:50,[8 14 44 45 47 48 49]);
time=time(off_ind);
ang=ang(off_ind);
data=data(:,off_ind,:);
data(isnan(data))=0;
data(isinf(data))=0;
[~,ang_ind]=sort(ang);
ang=ang(ang_ind);
data=data(:,ang_ind,:);
slice_tot = [7 14 15];
%%%================================
data_h=[];
ang_rate=1;
tau_rate=1;
thetan=ang;
for ele=1:size(data,1)
    data_h(ele,:,:)=data(ele,1:ang_rate:end,1:tau_rate:end);
    % for n=1:length(ang)
    %     data_h(ele,n,:)=map1D(data_h(ele,n,:),[0,1]);
    % end
end
thetan_real = thetan(1:ang_rate:end)'; 
thetan=mod(thetan_real,360);%linspace(-180,180,numThetan)+360;% must be positive.

if(ndims(data_h)==2)
    data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
end
% data_h=data_h(:,:,2:end);
numThetan=size(data_h,2);
nTau=size(data_h,3)-1;
N = nTau+1;%ceil((nTau+1)/sqrt(2)); 
%%==============================================
data_xrf_decom=[];
for ele=1:length(slice_tot)
data_xrf_decom(ele,:,:)=data_h(slice_tot(ele),:,:);
end
% save('tomopy_test.mat','data_xrf_decom','data_xrt','ang');
% return;
XRF_decom=permute(data_xrf_decom(:,:,:),[2 3 1]);
%%================================================
iR_num=zeros(N(1),N(1),NumElement);
