% ===============================================================================
% Tutorial for FAIR: Non-Parametric Image Registration
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% 
%   - data                 Hand, Omega=(0,20)x(0,25), 
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - transformation       affine2D
%   - regularizer          mbElastic
%   - optimizer            Gauss-Newton
% ===============================================================================

clear; close all; help(mfilename)

setupHandData
level = 5; omega = MLdata{level}.omega; m = MLdata{level}.m;
inter('reset','inter','splineInter','regularizer','moments','theta',0.01);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mbElastic','alpha',1000,'mu',1,'lambda',0);

[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
xc    = getStaggeredGrid(omega,m); % starting guess and reference for regularization
Rc    = inter(R,omega,center(xc,m)); 

% - initialize FAIR plots
FAIRplots('set','mode','NPIR-mb','fig',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% - run Non-Parametric Image Registration (Gauss-Newton)
fctn = @(yc) NPIRobjFctn(T,Rc,omega,m,xc,yc);  fctn([]);
yc = GaussNewton(fctn,xc,'Plots',@FAIRplots);
