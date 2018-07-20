Tomo_startup;
load PeriodicTable
filedir='./result/olga/';
% filedir='./';

files=dir([filedir,'opt3D_xy1_1.mat']); NumElement=2;

% files=dir('./result/miller/recon_shift_mc.mat'); NumElement=3;
Z=[26 34]; N=401;
% Z=[20 30 68]; N=64;
nslice=21;
n_delta=449*2;
for file=1:length(files)
    [pathstr,name,ext]=fileparts(files(file).name);
    output=[name,'.tiff'];
    load([filedir,name,ext]);
    recon=xCOR(n_delta+1:end);
    recon0=x(1:end);
    % recon=recon(1:N^2*NumElement*nslice,:);
    rr=reshape(full(recon),N,N,NumElement,nslice);
    rr0=reshape(full(recon0),N,N,NumElement,nslice);
    for ele=1:NumElement
        rr_temp=map1D(squeeze(rr(:,:,ele,:)),[0 1]);%;
        rr_temp0=map1D(squeeze(rr0(:,:,ele,:)),[0 1]);%;
        for k=1:nslice
            output=[name,'_',char(Element(Z(ele))),'.tiff'];
            output0=[name,'_0',char(Element(Z(ele))),'.tiff'];
            imwrite(rr_temp(:,:,k),output,'WriteMode','append');
            imwrite(rr_temp0(:,:,k),output0,'WriteMode','append');
        end
    end
end
