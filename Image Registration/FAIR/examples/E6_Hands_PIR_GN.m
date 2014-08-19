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
%   - pre-registration     rotation2D
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);

% setup data, interpolation, transformation, distance, 
% extract data from a particular level amd compute coefficients for interpolation
% setupHandData;
 viewImage('reset','viewImage','viewImage2D','colormap',gray(256),'axis','off');
% inter('reset','inter','splineInter','regularizer','moments','theta',1e-1);
% level =5; omega = MLdata{level}.omega; m = MLdata{level}.m;
% [T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega,'out',0);
% ===============================================================================
% T=double(imread('LenaRotate_P30.tiff'));
%T=double(imread('LenaCroppedRotate.tiff'));
% T=double(imread('LenaEye.tiff'));
% R=double(imread('LenaReference.tiff'));
T=double(imread('LoganRotate.tiff'));
R=double(imread('LoganReference.tiff'));
inter('reset','inter','splineInter','regularizer','moments','theta',1e-1);
omega=[0 size(R,1) 0 size(R,2)];
m = floor(size(R)/1);
T = inter('coefficients',T,[],omega,'regularizer','moments','theta',1e1);
R = inter('coefficients',R,[],omega,'regularizer','moments','theta',1e1);

% ===============================================================================
distance('reset','distance','SSD');
center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rigid2D','c',center); 
w0 = trafo('w0'); 
beta = 0; 
 M =[]; wRef = []; % disable regularization
%M = 2e2*speye(length(w0)); wRef = w0;     % initilize the regularization

% initialize plots
FAIRplots('reset','mode','PIR-GN','fig',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 
pause;
xc = getCellCenteredGrid(omega,m); 
Rc = inter(R,omega,xc);
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc);

% optimize
% w0=[0.1 200 -30]';
% w0=[0.4 260 -7]';
%w0=[0.1 90 0]';
 [wc,his] = GaussNewton(fctn,w0,'Plots',@FAIRplots,'solver','direct');
% [wc,his] = lBFGS(fctn,w0,'Plots',@FAIRplots,'solver','direct');