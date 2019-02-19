load PeriodicTable
slice=1;
sample='Paunesku';
setup_paunesku;
filedir=['./result/',sample,'/'];

files=dir([filedir,'Recon_shift_mc_weighted_TV.mat']); 
% files=dir([filedir,'opt3D_xy_weighted1_1FenoB.mat']); 
NumElement=6;
N_delta=numThetan;
for file=1:length(files)
    [pathstr,name,ext]=fileparts(files(file).name);
    output=[name,'.tiff'];
    load([filedir,name,ext]);
    recon=recon(N_delta+1:end,:);
    rr=reshape(full(recon),N,N,NumElement,nslice);
    % rr0=reshape(full(recon0),N,N,NumElement,nslice);
    for ele=1:NumElement
        rr_temp=map1D(squeeze(rr(:,:,ele,:)),[0 1]);%;
        % rr_temp0=map1D(squeeze(rr0(:,:,ele,:)),[0 1]);%;
        for k=1:nslice
            output=[name,'_',char(Element(Z(ele))),'.tiff'];
            % output0=[name,'_0',char(Element(Z(ele))),'.tiff'];
            imwrite(rr_temp(:,:,k),output,'WriteMode','append');
            % imwrite(rr_temp0(:,:,k),output0,'WriteMode','append');
        end
    end
end
