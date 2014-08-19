close all;
clear all;
Ref=imread('LoganReference.tiff');
omega=[0 size(Ref,1) 0 size(Ref,2)];
m=size(Ref);
x = getCellCenteredGrid(omega,m);
Crop=[50 20 50 50];
T=imcrop(flipud(Ref),Crop);
imshow(flipud(T))
return;
omegat=[Crop(1) Crop(1)+Crop(3)+1 Crop(2) Crop(2)+Crop(4)+1];
mt=size(T);
xc = getCellCenteredGrid(omegat,mt);
%%============ Affine Transformation
%%%------- Rigid
alpha=pi/6;
R = [cos(alpha),-sin(alpha);sin(alpha),cos(alpha)];
center = (omegat(2:2:end)+omegat(1:2:end))'/2;
w3=0;
w6=0;
shift=(eye(2)-R)*center;
trans=shift+[w3;w6];
trans=[w3;w6];
wc=[cos(alpha);-sin(alpha);trans(1);sin(alpha);cos(alpha);trans(2)];
%%%------------------------------------------------------------------
wc=[cos(alpha) -sin(alpha) shift(1) sin(alpha) cos(alpha) shift(2)]';
yc=affine2D(wc,xc);
%%%==============================
Tc = linearInter(T,omegat,yc);
yy=reshape(yc,length(yc)/2,2);
ytot=[reshape(yc,length(yc)/2,2);reshape(x,length(x)/2,2)];
yt=[reshape(yc,length(yc)/2,2);reshape(xc,length(xc)/2,2)];
Omega=[min(ytot(:,1)),max(ytot(:,1)),min(ytot(:,2)),max(ytot(:,2))];
   figure(1); 
%    viewImage2D(Ref,omega,m,'colormap','gray(256)');
%  hold on;
%  plotGrid(x,omega,m,'spacing',ceil(m/32));axis(omega);hold off
plotGrid(xc,omegat,mt,'spacing',floor(mt/5));
hold on;
 plotGrid(yc,omegat,mt,'spacing',floor(mt/5),'color','r');
axis(Omega);
 hold off
 figure(2); 
 viewImage2D(T,omegat,mt,'colormap','gray(256)');
 axis(omega)
 figure(3)
 viewImage2D(Tc,omegat,mt,'colormap','gray(256)');
 figure(4)
imshow(T,yy(:,1),yy(:,2))
axis xy equal
%  hold on;
%  plotGrid(yc,omegat,mt,'spacing',floor(mt/20),'color','r');
% 
% ytot1=[reshape(yc,length(yc)/2,2);reshape(xc,length(xc)/2,2)];
% 
% axis([min(ytot1(:,1)),max(ytot1(:,1)),min(ytot1(:,2)),max(ytot1(:,2))])
% %axis(omega);
% hold off;
% %title('w=[cos(pi/6);-sin(pi/6);0;sin(pi/6);cos(pi/6);0]','FontSize',13,'FontWeight','bold');
% title('w=[1; 0; 10.7; 0; 1; -22.5 ]','FontSize',18,'FontWeight','bold');