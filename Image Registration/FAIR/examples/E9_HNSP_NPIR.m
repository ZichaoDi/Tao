% ===============================================================================
% Example for NPIR, Non-Parametric Image Registration
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% 
%   - data                 HNSP, Omega=(0,2)x(0,1), level=4, m=[32,16]
%   - viewer               viewImage2D
%   - interpolation        linearInter
%   - distance             SSD
%   - regularizer          mbElastic
%   - optimizer            Gauss-Newton
% ===============================================================================

% setup data and initialize image viewer
setupHNSPData; FAIRdiary; level = 4; omega = MLdata{level}.omega; m = MLdata{level}.m; 

% initialize the interpolation scheme and coefficients
inter('reset','inter','splineInter'); 
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega,'out',0);
xc    = getCellCenteredGrid(omega,m); 
Rc    = inter(R,omega,xc);

% initialize distance measure
distance('set','distance','SSD');       

% initialize regularization, note: yc-yRef is regularized, elastic is staggered 
regularizer('reset','regularizer','mbElastic','alpha',1e4,'mu',1,'lambda',0);
y0   = getStaggeredGrid(omega,m); yRef = y0; yStop = y0;


% setup and initialize plots 
FAIRplots('reset','mode','NPIR-Gauss-Newton','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% build objective function, note: T coefficients of template, Rc sampled reference
fctn = @(yc) NPIRobjFctn(T,Rc,omega,m,yRef,yc); fctn([]); % report status

% -- solve the optimization problem -------------------------------------------
[yc,his] = GaussNewton(fctn,y0,'maxIter',500,'Plots',@FAIRplots,'yStop',yStop);
% report results
iter = size(his.his,1)-2; reduction = 100*fctn(yc)/fctn(y0);
fprintf('reduction = %s%% after %d iterations\n',num2str(reduction),iter);
diary off
