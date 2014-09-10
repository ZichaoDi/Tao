% function [R,X]=Radon3D;
close all;
% E=[100 0.5 0.5 0.5 0 0 0 0 0 0];
E =    [ ...
    0.08  .6900  .920  .810      0       0       0      0      0      0
%     -.8  .6624  .874  .780      0  -.0184       0      0      0      0
    0.1  .1100  .310  .220    .22       0       0    -18      0     10
    0.05  .1600  .410  .280   -.22       0       0     18      0     10
    0.09  .2100  .250  .410      0     .35    -.15      0      0      0
%     .1  .0460  .046  .050      0      .1     .25      0      0      0
%     .1  .0460  .046  .050      0     -.1     .25      0      0      0
%     .1  .0460  .023  .050   -.08   -.605       0      0      0      0
%     .1  .0230  .023  .020      0   -.606       0      0      0      0
%     .1  .0230  .046  .020    .06   -.605       0      0      0      0
    ];

I=phantom3d(E); % create the 3D phantom
theta=0:2:180;
[r1,X]=radon(squeeze(I(1,:,:)),theta);
R=zeros(size(X,1),size(I,1),length(theta));
R(:,1,:)=r1;
for i=2:size(I,1)
    [r]=radon(squeeze(I(i,:,:)),theta);
    R(:,i,:)=r;
end
figure('Name','Original 3D sample')
imshow3D(I)
colormap hot
figure('Name','2D projection (theta(i))'),
imshow3D(R) % view each 2D projections for different angle
colormap hot
figure('Name','sinogram(theta) v.s. slice'),
imshow3D(permute(R,[1 3 2])); % view each 1D projection for different slice
colormap hot
IP=iradon(squeeze(R(:,1,:)),theta);
Ip=zeros([size(I,1),size(IP)]);
for i=1:size(R,2)
    IP=iradon(squeeze(R(:,i,:)),theta);
    Ip(i,:,:)=reshape(IP,[size(IP),1]);
end
figure('Name','Reconstruced sample'),
imshow3D(Ip)
colormap hot