load('projs_large_fov.mat','projs_align_y');
load('projs_large_fov.mat','angles');
[nslice,ntau,ntheta]=size(projs_align_y);
sub_theta=[1:ntheta];
sub_tau=1:ntau;%[201:600];%  
% thetan=thetan(sub_theta);
[ang,ang_ind]=sort(angles(sub_theta));
whos projs_align_y
slice_sub=[1:10:1171];
nslice = length(slice_sub);
data=projs_align_y(slice_sub,sub_tau,sub_theta(ang_ind));
for i=1:nslice
    d=squeeze(data(i,:,:));
    d=d(:);
    num=[901,161];
    save(['test_chip_full',num2str(i),'.dat'],'num','-ascii');
    save(['test_chip_full',num2str(i),'.dat'],'d','-ascii','-append');
end
return;

%%%================================
NumElement=1;
ang_rate=1;
tau_rate=1;
data_h=data(:,1:tau_rate:end,1:ang_rate:end);
thetan =ang(1:ang_rate:end);% linspace(1,360,48);% 
clear projs_align_y angles data
thetan=thetan';
% save('angles.dat','thetan','-ascii');
whos data_h
return;
numThetan=length(thetan);
nTau=size(data_h,2)-1;
N = nTau+1;%ceil((nTau+1)/sqrt(2)); 

%%================================================
iR_num=zeros(N,N,NumElement);
nchannel=nslice;
Z=13;
