close all;
clear all;
Ref=imread('LenaReference.tiff');
omega=[0 size(Ref,1) 0 size(Ref,2)];
m=size(Ref);
x = getCellCenteredGrid(omega,m);
Crop=[200 200 100 100];
T=imcrop(Ref,Crop);
omegat=[Crop(1) Crop(1)+Crop(3)+1 Crop(2) Crop(2)+Crop(4)+1];
mt=size(T);
xc = getCellCenteredGrid(omega,m);
%%============ Affine Transformation
%%%------- Rigid
alpha=pi/6;
R = [cos(alpha),-sin(alpha);sin(alpha),cos(alpha)];
center = (omegat(2:2:end)+omegat(1:2:end))'/2;
w3=-49.5;
w6=50.5;
trans=(eye(2)-R)*center+[w3;w6];
wc=[cos(alpha);-sin(alpha);trans(1);sin(alpha);cos(alpha);trans(2)];
%%%------------------------------------------------------------------
%  wc=[1 -0.2 25.5 0.2 0.75 0]';
yc=affine2D(wc,xc);
%%%==============================
Tc = linearInter(Ref,omega,yc);
   figure(1); 
%   viewImage2D(Ref,omega,m,'colormap','gray(256)');
%  hold on;
%  plotGrid(x,omega,m,'spacing',ceil(m/32));axis(omega);hold off
plotGrid(xc,omega,m,'spacing',floor(m/20));
hold on;
 plotGrid(yc,omega,m,'spacing',floor(m/20),'color','r');
 ytot=[reshape(yc,length(yc)/2,2);reshape(xc,length(xc)/2,2)];
axis([min(ytot(:,1)),max(ytot(:,1)),min(ytot(:,2)),max(ytot(:,2))])

 hold off
 figure(2); 
 viewImage2D(T,omegat,mt,'colormap','gray(256)');
 axis(omega)
 figure(3)
 viewImage2D(Tc,omega,m,'colormap','gray(256)');
 hold on;
 plotGrid(yc,omega,m,'spacing',floor(m/20),'color','r');
axis([min(ytot(:,1)),max(ytot(:,1)),min(ytot(:,2)),max(ytot(:,2))])
% %axis(omega);
% hold off;
% %title('w=[cos(pi/6);-sin(pi/6);0;sin(pi/6);cos(pi/6);0]','FontSize',13,'FontWeight','bold');
% title('w=[1; 0; 10.7; 0; 1; -22.5 ]','FontSize',18,'FontWeight','bold');