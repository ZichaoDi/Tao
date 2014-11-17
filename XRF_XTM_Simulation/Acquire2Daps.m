XTM=[];
Angle=csvread('glass_rod_W_Au_wire_Scan_0.csv',1,1);
p1=37; 
thetan=Angle(p1:end);
for i=p1-1:length(Angle)-1
data=h5read(['/2xfm0415/2xfm_0',num2str(i+100),'.h50'],'/MAPS/scalers');
XTM_ic=data(:,:,3);
angle=Angle(i);
x=data(:,:,17);
y=data(:,:,18);
XTM=XTM_ic;
XRF=h5read(['/2xfm0415/2xfm_0',num2str(i+100),'.h50'],'/MAPS/mca_arr');
DetChannel=h5read(['/2xfm0415/2xfm_0',num2str(i+100),'.h50'],'/MAPS/energy');
elements=h5read(['/2xfm0415/2xfm_0',num2str(i+100),'.h50'],'/MAPS/channel_names');
save(['data/2xfm0415/2xfm_0',num2str(i+100),'.mat'],'x','y','elements','DetChannel','XTM','XRF','angle');
end



xx=cell(20,1);
for i=1:20
load(['2xfm_0',num2str(i+135),'.mat']);
xx{i}=XTM;
end