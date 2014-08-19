%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: 1D interpolation, scale
%
% the tutorial explores the scale-space idea for 1D data
% the spline based model is based on the minimization of
% D(T)+theta S(T)!= min, 
% where D(T) is a data fitting term and S(T) is the linearized bending energy
%
% - setup 1D data
% - visualize data
% - compute approximations (spline, various theta's) and visualize
%==============================================================================

clear, close all, help(mfilename); echo on

%% setup data 
fprintf('%s\n','generate noisy data ');

omega = [0,10]; m = 21;
dataX = getCellCenteredGrid(omega,m);
dataT  = rand(m,1); dataT([1,end]) = 0; % T should be compactly supported

% visualize the data
figure(2); clf;
ph = plot(dataX,0*dataT,'k.',dataX,dataT,'b.','markersize',30); hold on; 
title('1D multiscale','fontsize',20);
axis([omega(1)-1,omega(2)+1,min(dataT)-1,max(dataT)+1]);
lstr = {'location','data'}; legend(ph,lstr,1);
pause;

%% discretization for the model  -----------------------------------------------
m  = 101; % discretization for the model, h=omega./m
xc = getCellCenteredGrid(omega,m);

% the coefficients depend on the scale-parameter theta
T = @(theta) getSplineCoefficients(dataT,...
  'regularizer','moments','theta',theta,'dim',1,'out',0);

% the continuous model
Tc = @(theta) splineInter(T(theta),omega,xc);

% show results for various theta's: fine to coarse scale
fprintf('%s\n','show results for various theta''s');
theta = [0,logspace(-3,3,51)];
for j=1:length(theta),
  fprintf('.');
  if j>1, 
    set(qh,'visible','off');  
  end;
  qh = plot(xc,Tc(theta(j)),'g-','linewidth',3);
  ylabel(sprintf('\\theta=%s',num2str(theta(j))),'fontsize',20); pause(1/10)
  if j==1,  lstr{3} = 'spline';  legend([ph;qh],lstr,1); pause; end;
end;
fprintf('\n'); pause;

% show results for various theta's: coarse to fine scale
for j=length(theta):-1:1,
  fprintf('.');
  set(qh,'visible','off');
  qh = plot(xc,Tc(theta(j)),'g-','linewidth',3);
  ylabel(sprintf('\\theta=%s',num2str(theta(j))),'fontsize',20); pause(1/10)
  if j==1, pause; end;
end;
fprintf('\n'); 

echo off