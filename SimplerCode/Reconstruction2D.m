function Reconstruction2D;
close all
% I=ones(40,40);
% I(2,2)=0.5;
% I(2,3)=0.7;
% I(3,2)=1;
% I(3,3)=0.9;

x=linspace(0,100,100);
y=linspace(0,100,100);
[X,Y] = meshgrid(x,y);
I=(X-50).^2+(Y-50).^2;
figure('name','Object');
imagesc(I)
theta=45;%[0 45];%0:180;%
[R,xp]=radon(I,theta);
R1=R(:,1);
figure('name','Projection')
plot(1:length(R1),R1,'r.-')
figure('name','original sample')
imagesc(I);
figure('name','Radon Transform')
imagesc(theta,xp,R);
Ip=iradon(R,theta);
figure('name','Reconstructed sample');
imagesc(Ip)
TransR=exp(-R);
