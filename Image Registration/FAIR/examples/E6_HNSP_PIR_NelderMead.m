% ===============================================================================
% (c) Jan Modersitzki 2010-08-11, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki/software/fair.html}
%
% Tutorial for FAIR: PIR,  Parametric Image Registration
% 
%   - data                 Hand, Omega=(0,20)x(0,25), level=5, m=[32,32]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - optimization         Nelder-Mead (fminsearch)
% ===============================================================================

clear, close all, help(mfilename);

% setup data, interpolation, transformation, distance, 
% extract data from a particular level amd compute coefficients for interpolation
setupHNSPData; 
inter('set','inter','splineInter'); 
level = 4; omega = MLdata{level}.omega; m = MLdata{level}.m; 
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
xc = getCellCenteredGrid(omega,m); 
Rc = inter(R,omega,xc);
distance('set','distance','SSD');       
trafo('reset','trafo','affine2D');
w0 = trafo('w0'); 

% build objective function
% note: T  is template image
%       Rc is sampled reference
%       optional Tikhonov-regularization is disabled by setting m = [], wRef = []
%       beta = 0 disables regularization of Hessian approximation
beta = 0; M = []; wRef = [];
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc); 
fctn([]);   % report status

optn = optimset('display','iter');
[ycNM,JcNM,exitFlag,out] = fminsearch(fctn,w0,optn);
yc = trafo(ycNM,xc);
showResults(MLdata,yc,'level',level)
