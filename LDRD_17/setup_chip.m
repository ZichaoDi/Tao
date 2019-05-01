ang_rate=10;
tau_rate=2;
nslice = 870;
load('data/projs_small_fov.mat','projs_align_y');
load('data/projs_small_fov.mat','angles');
[~,ang_ind]=sort(angles);
ang=angles(ang_ind);
slice_sub=[100];
data=projs_align_y(slice_sub,:,ang_ind);
NumElement=1;
%%%================================
thetan=ang;
data_h=data(:,2:tau_rate:end,1:ang_rate:end);
clear projs_align_y angles data

thetan_real =thetan(1:ang_rate:end)';% linspace(1,360,48);% 
thetan=mod(thetan_real,360);%linspace(-180,180,numThetan)+360;% must be positive.
numThetan=length(thetan);
nTau=size(data_h,2)-1;
N = nTau+1;%ceil((nTau+1)/sqrt(2)); 

%%================================================
iR_num=zeros(N,N,NumElement);
