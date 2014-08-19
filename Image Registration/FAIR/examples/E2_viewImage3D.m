%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: load and visualize 3D MRI data
%
% - load data             ('US.jpg')
% - setup viewer          (viewImage2D)
% - view  data
%==============================================================================

clear, close all, help(mfilename), 
echo on


% setup data -----------------------------------------------------------------
fprintf('%s\n','generate data (MRI)');
load mri;
omega = [0,20,0,20,0,10];
T     = double(squeeze(D));  % specify physical domain
m     = size(T);

%%
fprintf('%s\n','use image montage for visualization')
FAIRfigure(1); clf; 
imgmontage(T,omega,m,'colormap',gray(100));
title('visualization of 3D data - imgmontage','fontsize',30)
pause;

%%
fprintf('%s\n','use a volumetric view')
FAIRfigure(2); clf; 
volView(T,omega,m,'facecolor',[240,140,100]/256,'facealpha',0.75); hold on;
colormap(gray(100))
pause;

%%
fprintf('%s\n','use a volumetric view + slices')
vh = viewSlices(T,omega,m);
pause;

figure(3); clf; colordef(gcf,'black');
volView(T,omega,m,'facecolor',[240,140,100]/256,'facealpha',0.75); hold on;
colormap(gray(100))
vh = viewSlices(T,omega,m);
ah = viewSlices(T,omega,m,'s1',64,'s2',64,'s3',[]);
pause; 

for j=1:m(3);
  set(vh,'visible','off');
  vh = viewSlices(T,omega,m,'s1',[],'s2',[],'s3',[j]);
  pause(1/10);
end;
echo off