% ===============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: PIR,  Parametric Image Registration
% 
%   - data                 Hand, Omega=(0,20)x(0,25), level=5, m=[32,32]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     rotation2D, not regularized
%   - optimization         Steepest Descent
% ===============================================================================

clear, close all, help(mfilename);

% setup data, interpolation, transformation, distance, 
% extract data from a particular level amd compute coefficients for interpolation
setupHandData;
inter('reset','inter','splineInter','regularizer','moments','theta',1e-2);
level = 5; omega = MLdata{level}.omega; m = MLdata{level}.m;
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega,'out',0);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D'); 
w0 = trafo('w0'); beta = 0; M =[]; wRef = []; % disable regularization

% initialize plots
FAIRplots('reset','mode','PIR-SD','fig',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 


xc = getCellCenteredGrid(omega,m); 
Rc = inter(R,omega,xc);
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc);

% run STEEPEST DESCENT with an initial stepLength
[Jc,para,dJ] = fctn(w0); stepLength = 0.5/norm(dJ);
optn = {'stepLength',stepLength,'maxIter',50,'Plots',@FAIRplots};
wc =  SteepestDescent(fctn,w0,optn{:});

