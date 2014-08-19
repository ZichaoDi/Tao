% ===============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% (Counter) Example for EPI registration. Please try to register
% both images that show deformations into opposite directions. 
% 
%   - data                 2D EPI slices m=[256,256]
%   - viewer               viewImage2D
%   - interpolation        splineInterMex
%   - distance             SSD
%   - optimization         Gauss-Newton
%
% see also apps/EPI/contents
% ===============================================================================

clear; clc; %close all
setup2DEPIData;

level = 5; omega = MLdata{level}.omega; m = MLdata{level}.m; 
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
xc = getCellCenteredGrid(omega,m); 
Rc = inter(R,omega,xc);


% initialize the transformation and a starting guess
trafo('reset','trafo','splineTransformation2Dsparse','omega',omega,'p',4*[3 4],'m',m);
%trafo('reset','trafo','affine2Dsparse','omega',omega,'p',[3 4],'m',m);
w0 = trafo('w0');


% setup plots and initialize
FAIRplots('reset','mode','PIR-GN','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% build objective function
% note: T  is data for template image
%       Rc is sampled reference image
%       optional Tikhonov-regularization is disabled by setting m = [], wRef = []
%       beta = 0 disables additiona regularization of Hessian approximation
beta = 1; M = 1e5*speye(numel(w0)); wRef = 0;
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc); 
fctn([]);   % report status
[a,para,da] = fctn(w0)

%% -- solve the optimization problem on one level
% [wc,his] = GaussNewton(fctn,w0,'Plots',@FAIRplots); return;


%% finally: run the MultiLevel Non-Parametric Image Registration
[wc,his] = MLPIR(MLdata,'minLevel',6,'M',M,'beta',beta,'wRef',wRef);