% Tutorial for FAIR: Sparse spline transformation 3D
% 
%           data: 3D cardiac gated mice data
% transformation: splineTransformation3Dsparse

clear; clc; %close all
setup3DmiceData;
inter('reset','inter','linearInter');


level = 6; omega = MLdata{level}.omega; m = MLdata{level}.m; 
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
xc = getCellCenteredGrid(omega,m); 
Rc = inter(R,omega,xc);


% initialize the transformation and a starting guess
% trafo('reset','trafo','affine3Dsparse');
 trafo('reset','trafo','splineTransformation3Dsparse','omega',omega,'p',[3 4 3],'m',m);
%trafo('reset','trafo','affine3Dsparse','omega',omega,'p',[3 4 3],'m',m);
w0 = trafo('w0');


% setup plots and initialize
FAIRplots('reset','mode','PIR-GN','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% build objective function
% note: T  is data for template image
%       Rc is sampled reference image
%       optional Tikhonov-regularization is disabled by setting m = [], wRef = []
%       beta = 0 disables additiona regularization of Hessian approximation
beta = 0; M = []; wRef = [];
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc); 
fctn([]);   % report status
[a,para,da] = fctn(w0)

%% -- solve the optimization problem on one level
dummy = @(a,b) []; % @FAIRplots
% [wc,his] = GaussNewton(fctn,w0,'Plots',dummy,'maxIter',5); return;


%% finally: run the MultiLevel Non-Parametric Image Registration
tic;
[wc,his] = MLPIR(MLdata,'minLevel',4);

toc;
