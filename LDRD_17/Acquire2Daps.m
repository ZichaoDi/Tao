% XTM=[];
% Angle=csvread('glass_rod_W_Au_wire_Scan_0.csv',1,1);
IndPro=173:245;
thetan=zeros(length(IndPro),1);
for i=1:length(IndPro);%p1-1:length(Angle)-1
data=h5read(['/nfs2/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/ApsDataExtract/2xfm_SV/2xfm_0',num2str(IndPro(i)),'.h5'],'/MAPS/scalers');
XTM_ic=data(:,:,4);
% angle=str2num(aa(27:end));
% thetan(i)=angle;
x=data(:,:,17);
y=data(:,:,18);
XTM=XTM_ic;
XRF=h5read(['/nfs2/wendydi/Documents/Research/Tao/XRF_XTM_Simulation/data/ApsDataExtract/2xfm_SV/2xfm_0',num2str(IndPro(i)),'.h5'],'/MAPS/mca_arr');
% XRF=h5read(['/Users/Wendydi/Documents/MATLAB/APSdata/2xfm1211_14/2xfm_0',num2str(IndPro(i)),'.h5'],'/MAPS/mca_arr');
% DetChannel=h5read(['/Users/Wendydi/Documents/MATLAB/APSdata/2xfm1211_14/2xfm_0',num2str(IndPro(i)),'.h5'],'/MAPS/energy');
% elements=h5read(['/Users/Wendydi/Documents/MATLAB/APSdata/2xfm1211_14/2xfm_0',num2str(IndPro(i)),'.h5'],'/MAPS/channel_names');
% save(['data/2xfm1211_14_4/2xfm_',num2str(IndPro(i)),'.mat'],'x','y','elements','DetChannel','XTM','XRF','angle');
save(['data/2xfm1211_14_4/2xfm_',num2str(IndPro(i)),'.mat'],'XTM');

end



