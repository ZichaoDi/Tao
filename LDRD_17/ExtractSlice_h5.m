slice=30;
nTau=1750;
numThetan=73;
numMAPSchannel=34;
data=zeros(numMAPSchannel,numThetan, nTau);
for t=1:numThetan
    data_t=h5read(['/nfs2/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/ApsDataExtract/2xfm_SV/2xfm_0',num2str(172+t),'.h5'],'/exchange_4/data');
    data(:,t,:)=permute(data_t(:,slice,:),[3 2 1]);
end
channel_names=h5read(['/nfs2/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/ApsDataExtract/2xfm_SV/2xfm_0173.h5'],'/exchange_4/channel_names');
save(['data/ApsDataExtract/2xfm_SV/slice',num2str(slice),'.mat'],'data','channel_names');
