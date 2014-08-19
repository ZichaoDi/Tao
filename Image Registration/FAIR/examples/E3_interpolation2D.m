%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: interpolation in 2D
%
% - load data   (setupUSData)
% - interpolate (linearInter) on different resolutions
% - visualize   (viewImage2D)
% - interpolate (splineInter) on different resolutions
% - visualize   (viewImage2D)
%==============================================================================

setupUSData; close all; 
T = dataT; xc = @(m) getCellCenteredGrid(omega,m); 

inter('set','inter','linearInter');
for p=5:7,
  m = 2^p*[1,1]; 
  Tc = inter(T,omega,xc(m));
  figure(p-4); viewImage2D(Tc,omega,m); colormap(gray(256));
end;

inter('set','inter','splineInter');
T = getSplineCoefficients(dataT,'regularizer','moments','theta',100);
for p=5:7,
  m = 2^p*[1,1]; 
  Tc = inter(T,omega,xc(m));
  figure(p-1); viewImage2D(Tc,omega,m); colormap(gray(256));
end;
