% ===============================================================================
% Example for Non-Parametric Image Registration with pre-registration
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% 
%   - data                 HNSP, Omega=(0,2)x(0,1), m=[ 512 256], level = 4
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - regularizer          mfElastic
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


% -- the PIR pre-registration -------------------------
% intialize the pre-registration
trafo('reset','trafo','rigid2D'); w0 = trafo('w0')
% setup plots and initialize objective function for PIR
FAIRplots('set','mode','PIR-GN-rigid','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 
beta = 0; M = []; wRef = []; xc = getCellCenteredGrid(omega,m);
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc); fctn([]); % report status
% solve the PIR
[wc,his] = GaussNewton(fctn,w0,'maxIter',500,'Plots',@FAIRplots);

% prolongate intermediates to level = 7
level = 5; omega = MLdata{level}.omega; m = MLdata{level}.m; 
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
xc    = getCellCenteredGrid(omega,m); 
Rc    = inter(R,omega,xc);
yRef  = getStaggeredGrid(omega,m); yStop = yRef;
y0    = grid2grid(trafo(wc,getNodalGrid(omega,m)),m,'nodal','staggered');
% -- END the PIR pre-registration ---------------------

% setup and initialize plots 
FAIRplots('set','mode','NPIR-Gauss-Newton','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% build objective function, note: T coefficients of template, Rc sampled reference
fctn = @(yc) NPIRobjFctn(T,Rc,omega,m,yRef,yc); fctn([]); % report status

% -- solve the optimization problem -------------------------------------------
[yc,his] = GaussNewton(fctn,y0,'maxIter',500,'Plots',@FAIRplots,'yStop',yStop);
% report results
iter = size(his.his,1)-2; reduction = 100*fctn(yc)/fctn(y0);
fprintf('reduction = %s%% after %d iterations\n',num2str(reduction),iter);
diary off
