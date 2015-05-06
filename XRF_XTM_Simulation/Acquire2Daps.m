% XTM=[];
% Angle=csvread('glass_rod_W_Au_wire_Scan_0.csv',1,1);
IndPro=173:245;
thetan=zeros(length(IndPro),1);
for i=1:length(IndPro);%p1-1:length(Angle)-1
% data=h5read(['/Users/Wendydi/Documents/MATLAB/APSdata/2xfm1211_14/2xfm_0',num2str(IndPro(i)),'.h5'],'/MAPS/scalers');
% % temp=h5read(['/Users/Wendydi/Documents/MATLAB/APSdata/2xfm1211_14/2xfm_0',num2str(IndPro(i)),'.h5'],'/MAPS/extra_strings');
% % aa=temp{end};
% XTM_ic=data(:,:,4);
% % angle=str2num(aa(27:end));
% % thetan(i)=angle;
% x=data(:,:,17);
% y=data(:,:,18);
% XTM=XTM_ic;
% XRF=h5read(['/Users/Wendydi/Documents/MATLAB/APSdata/2xfm1211_14/2xfm_0',num2str(IndPro(i)),'.h5'],'/MAPS/mca_arr');
% DetChannel=h5read(['/Users/Wendydi/Documents/MATLAB/APSdata/2xfm1211_14/2xfm_0',num2str(IndPro(i)),'.h5'],'/MAPS/energy');
% elements=h5read(['/Users/Wendydi/Documents/MATLAB/APSdata/2xfm1211_14/2xfm_0',num2str(IndPro(i)),'.h5'],'/MAPS/channel_names');
% save(['data/2xfm1211_14_4/2xfm_',num2str(IndPro(i)),'.mat'],'x','y','elements','DetChannel','XTM','XRF','angle');
XRF_roi=h5read(['/Users/Wendydi/Documents/MATLAB/APSdata/2xfm1211_14/2xfm_0',num2str(IndPro(i)),'.h5'],'/MAPS/XRF_roi');

save(['data/2xfm1211_14/2xfm_',num2str(IndPro(i)),'.mat'],'XRF_roi','-append');

% save(['data/2xfm1211_14_4/2xfm_',num2str(IndPro(i)),'.mat'],'XTM');

end



% xx=cell(20,1);
% for i=1:20
% load(['2xfm_0',num2str(i+135),'.mat']);
% xx{i}=XTM;
% end