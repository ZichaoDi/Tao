Tomo_startup;
load PeriodicTable
filedir='./result/paunesku1/';

files=dir([filedir,'Recon_shift_weighted.mat']); NumElement=6;

% files=dir('./result/miller/recon_shift_mc.mat'); NumElement=3;
Z=[15 16 20 22 26 30]; N=81;
% Z=[20 30 68]; N=64;
nslice=131;
for file=1:length(files)
    [pathstr,name,ext]=fileparts(files(file).name);
    output=[name,'.tiff'];
    load([filedir,name,ext]);
    recon=recon(end:-1:1,:);
    recon=recon(1:N^2*NumElement,:);
    recon=recon(end:-1:1,:);
    rr=reshape(full(recon),N,N,NumElement,nslice);
    for ele=1:NumElement
        rr_temp=squeeze(rr(:,:,ele,:));
        for k=1:nslice
            output=[name,'_',char(Element(Z(ele))),'.tiff'];
            % output=[name,'.tiff'];
            imwrite(uint16(rr_temp(:,:,k)),output,'WriteMode','append');
        end
    end
end
