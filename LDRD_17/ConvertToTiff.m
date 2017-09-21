
files=dir('./result/paunesku/*.mat');
N=83;
nslice=131;
for file=1:length(files)
    [pathstr,name,ext]=fileparts(files(file).name);
    output=[name,'.tiff'];
    load([name,ext]);
    recon=recon(end:-1:1,:);
    recon=recon(1:N^2,:);
    recon=recon(end:-1:1,:);
    rr=reshape(full(recon),N,N,nslice);
    for k=1:nslice
        imwrite(uint16(rr(:,:,k)),output,'WriteMode','append');
    end
end
