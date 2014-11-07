XTM=[];
Angle=csvread('glass_rod_W_Au_wire_Scan_0.csv',1,1);
p1=37; 
thetan=Angle(p1:end);
for i=p1-1:length(Angle)-1
data=h5read(['2xfm_0',num2str(i+100),'.h50'],'/MAPS/scalers');
XTM=[XTM; data(1000,:,4)];
end
XTM=XTM(:,:);
XTM(find(isnan(XTM)))=0;

thetan=thetan(1:end);
%  imagesc(XTM)