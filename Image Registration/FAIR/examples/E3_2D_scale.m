%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: 2D interpolation and visualization
%
% the tutorial explores the scale-space idea for 2D data
% the spline based model is based on the minimization of
% D(T)+theta S(T)!= min, 
% where D(T) is a data fitting term and S(T) is the linearized bending energy
%
% - load data ('US.jpg')
% - display and visualize data  (viewImage2D)
% - compute approximations (spline, various theta's) and visualize
%==============================================================================

clear, close all, help(mfilename); echo on

%% load data, define a doman and an initial discretization
dataT = double(imread('US.jpg'));
omega = [0,size(dataT,1),0,size(dataT,2)];
m     = [128,128]/2;
xc    = getCellCenteredGrid(omega,m);

% setup image viewer
viewImage('reset','viewImage','viewImage2D','colormap','gray(256)');

% setup spline interpolation
inter('reset','inter','splineInter','regularizer','moments');
inter('disp');

scaleParameter = [logspace(3,-2,23),0];
titleStr = @(j) title(...
  sprintf('scale-space, \\theta=%s',num2str(scaleParameter(j))),'fontsize',30);

for j=1:length(scaleParameter);
  % set scale parameter and compute spline coefficients
  inter('set','theta',scaleParameter(j));
  T = inter('coefficients',dataT,[],omega,'out',0);
  
  % interpolate image
  Tc = inter(T,omega,xc);
  
  if j==1, % initilize figure
    figure(1); clf;
    vh = viewImage(Tc,omega,m); titleStr(j);
    pause;
  else
    set(vh,'cdata',reshape(Tc,m)'); titleStr(j);
    pause(1/500);
  end;
end;
  
echo off