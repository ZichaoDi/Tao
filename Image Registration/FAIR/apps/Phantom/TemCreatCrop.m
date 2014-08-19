close all;
clear all;
warning off;
Ref=double(imread('LoganReference.tiff'));
Ref=flipud(Ref);
crop=[50 130 50 50];
[Rcrop,rect]=imcrop(uint8(Ref),crop);

%Rcrop=fliplr(flipud(Rcrop));
Rcrop=flipud(Rcrop);
omega=[-(size(Ref,1))/2 (size(Ref,1))/2 -(size(Ref,2))/2 (size(Ref,2))/2];
omegat=[crop(1) crop(1)+crop(3) crop(2) crop(2)+crop(4)]+[-(size(Ref,1))/2 -(size(Ref,1))/2 -(size(Ref,2))/2 -(size(Ref,1))/2];

figure
imshow(uint8(Rcrop))
figure
imshow(uint8(Ref))
%    return;
m=size(Ref);
xc = getCellCenteredGrid(omega,m);
viewImage2D(Rcrop,omegat,size(Rcrop),'colormap','gray(256)');
%   axis(omega)
%   hold on;
% plotGrid(xc,omega,m,'spacing',floor(m/10));
%  return;
%%%------- Rigid
alpha=-pi/6;
w1=cos(alpha);
w2=-sin(alpha);
w4=sin(alpha);
w5=cos(alpha);
w3=100;
w6=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%w1=1.5;w2=0;w3=0;w4=0.5;w5=1;w6=0;  %affine
%w1=1;w2=0;w3=130.5;w4=0;w5=1;w6=-90;  %Translation
tform=[w1 w4; w2 w5; w3 w6];
tform = maketform('affine',tform);
inv_points = tformfwd(tform, reshape(xc,prod(m),2));
bound=tformfwd(tform,[omega(1),omega(3);omega(2),omega(3);omega(1),omega(4);omega(2),omega(4)]);
sub_bound=tformfwd(tform,[omegat(1),omegat(3);omegat(2),omegat(3);omegat(1),omegat(4);omegat(2),omegat(4)]);
xx=reshape(xc,prod(m),2);
yy=inv_points;
OG=[min(yy(:,1)),max(yy(:,1)),min(yy(:,2)),max(yy(:,2))];
%OGx=[min(xx(:,1)),max(xx(:,1)),min(yt(:,2)),max(yy(:,2))];
yt=[yy;xx];
OGx=[min(yt(:,1)),max(yt(:,1)),min(yt(:,2)),max(yt(:,2))];

[Rc,xdata,ydata] = imtransform(Ref,tform,'UData',omega(1:2),'VData',omega(3:4),'XData',OG(1:2),'YData',OG(3:4),'XYScale',1); %,'UData',omega(1:2),'VData',omega(3:4),'XData',OG(1:2),'YData',OG(3:4)
figure, imshow(uint8(Ref))
figure,imshow(uint8(Rc))
figure,
viewImage2D(Ref,omega,m,'colormap','gray(256)');
hold on;
plot([omegat(1),omegat(2),omegat(1),omegat(2)],[omegat(3),omegat(3),omegat(4),omegat(4)],'bo')
figure,
viewImage2D(Rc,[xdata,ydata],size(Rc),'colormap','gray(256)');
hold on;
plot(sub_bound(:,1),sub_bound(:,2),'r*')
plot([xdata,xdata],[ydata(1),ydata(1),ydata(2),ydata(2)],'r.')

figure,
viewImage2D(Rc,OG,size(Rc),'colormap','gray(256)');
figure,
viewImage2D(Ref,omega,size(Ref),'colormap','gray(256)');
hold on;
plotGrid(xc,omega,m,'spacing',floor(m/10));
hold on;
plotGrid(yy(:),OGx,m,'spacing',floor(m/10),'color','r');
axis(OGx);
hold off
figure,
plotGrid(xc,omega,m,'spacing',floor(m/20));
hold on;
plotGrid(yy(:),omega,m,'spacing',floor(m/20),'color','r');
hold on;
% plot(sub_bound(:,1),sub_bound(:,2),'r*','MarkerSize',10)
% hold on;
% plot([omegat(1),omegat(2),omegat(1),omegat(2)],[omegat(1),omegat(3),omegat(4),omegat(4)],'b*','MarkerSize',10)
plot(sub_bound(1,1),sub_bound(1,2),'r*',sub_bound(2,1),sub_bound(2,2),'ro',sub_bound(3,1),sub_bound(3,2),'rs',sub_bound(4,1),sub_bound(4,2),'rd','MarkerFaceColor','r','MarkerSize',10);
hold on;
plot(omegat(1),omegat(3),'b*',omegat(2),omegat(3),'bo',omegat(1),omegat(4),'bs',omegat(2),omegat(4),'bd','MarkerFaceColor','b','MarkerSize',10);
hold off;
axis image
crop1=[min(sub_bound(:,1)),min(sub_bound(:,2)),max(sub_bound(:,1))-min(sub_bound(:,1)),max(sub_bound(:,2))-min(sub_bound(:,2))];
crop1(1:2)=crop1(1:2)-[xdata(1),ydata(1)];%-[w3,w6];
% sub_bound=sub_bound+repmat([size(Ref,1)/2,size(Ref,2)/2],4,1);
% xdata=xdata+size(Ref,1)/2;
% ydata=ydata+size(Ref,2)/2;
% figure,
% viewImage2D(Rc,[xdata,ydata],size(Rc),'colormap','gray(256)');
% hold on;
% plot(sub_bound(:,1),sub_bound(:,2),'r*')
% figure,
% imshow(uint8(Rc));
Rcrop=imcrop(uint8(Rc),crop1);
% imwrite(Rcrop,'LoganTest1.tiff');
figure,imshow(uint8(flipud(Rcrop)))
