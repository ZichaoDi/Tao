% ===============================================================================
% Example for NPIR, Non-oparametric Image Registration
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
% 
%   - data                 MRI (head), level=7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             MI
%   - pre-registration     affine2D
%   - regularizer          mbElastic
%   - optimization         lBFGS
% ===============================================================================

close all, help(mfilename);

setupMRIdata
inter('reset','inter','splineInter','regularizer','moments','theta',1e-2);
distance('reset','distance','MI','nT',8,'nR',8);
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mbElastic','alpha',1e-3,'mu',1,'lambda',0);


FAIRdiary; level = 7; omega = MLdata{level}.omega; m = MLdata{level}.m; 

% initialize the interpolation scheme and coefficients
inter('reset','inter','splineInter'); 
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega,'out',0);
xc    = getCellCenteredGrid(omega,m); 
Rc    = inter(R,omega,xc);

% initialize distance measure
distance('set','distance','MI','nT',8,'nR',8);       

% initialize regularization, note: yc-yRef is regularized, elastic is staggered 
regularizer('reset','regularizer','mbElastic','alpha',1e-2,'mu',1,'lambda',0);
y0   = getStaggeredGrid(omega,m); yRef = y0; yStop = y0;


% setup and initialize plots 
FAIRplots('set','mode','lBFGS','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 


% build objective function, note: T coefficients of template, Rc sampled reference
fctn = @(yc) NPIRBFGSobjFctn(T,Rc,omega,m,yRef,yc); fctn([]); % report status


% -- solve the optimization problem -------------------------------------------
[yc,his] = lBFGS(fctn,y0,'maxIter',500,'Plots',@FAIRplots,'yStop',yStop,'tolJ',1e-4);
% report results
% iter = size(his.his,1)-2; reduction = 100*fctn(yc)/fctn(getStaggeredGrid(omega,m)); yOpt = yc;
% fprintf('reduction = %s%% after %d iterations\n',num2str(reduction),iter);
diary off
% 
% [yc,wc,his] = MLIR(MLdata,'minLevel',4,'maxLevel',4,'parametric',1,'plotMLiter',0,'plotIter',1);
