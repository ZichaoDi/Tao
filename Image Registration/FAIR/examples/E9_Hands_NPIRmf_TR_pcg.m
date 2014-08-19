% Tutorial for FAIR: Non-Parametric Image Registration
% (c) Jan Modersitzki 2009/04/24, see FAIR.2 and FAIRcopyright.m.
% - initialize
%     data:            Hands,           see setupHandData, level = 5;
%     visualization:   viewImage2D,     see viewImage.m
%     interpolation:   splineInter,   see inter.m
%     distance:        SSD,             see distance.m
%     trafo:           affine2D,        see trafo.m
%     regularizer:     mfElastic,       see regularizer.m
% - initialize FAIR plots
% - run optimization
%     NPIR:            Non-Parametric Image Registration, Trust-Region, PCG

clear; close all; help(mfilename)

setupHandData
level = 5; omega = MLdata{level}.omega; m = MLdata{level}.m;
inter('reset','inter','splineInter','regularizer','moments','theta',0.01);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mfElastic','alpha',1000,'mu',1,'lambda',0);

[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
xc    = getStaggeredGrid(omega,m); % starting guess and reference for regularization
Rc    = inter(R,omega,center(xc,m)); 

% - initialize FAIR plots
FAIRplots('set','mode','TR-PCG-mf','fig',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% - run Non-Parametric Image Registration (Gauss-Newton)
fctn = @(yc) NPIRobjFctn(T,Rc,omega,m,xc,yc);  fctn([]);
yc = TrustRegion(fctn,xc,'pfun',@(x,H)MGsolver(x,H),'Plots',@FAIRplots);
