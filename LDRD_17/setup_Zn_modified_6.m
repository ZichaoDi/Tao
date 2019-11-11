data=[];
numThetan=22;
files=dir('~/Documents/Research/APS/Zn_modified_6/*.tif');
for n=1:numThetan
    [pathstr,name,ext]=fileparts(['~/Documents/Research/APS/Zn_modified_6/',files(n).name]);
    data(n,:,:)=imread([pathstr,'/', name, ext]);
end
XRF_decom3D=data;%permute(data,[1 3 2]);
nslice=size(XRF_decom3D,2);
nTau=size(XRF_decom3D,3)-1;
N=nTau+1;
thetan_real=importdata('angle.txt');
thetan_real=thetan_real';
Z=30;
