% ===============================================================================
% Example for MLSIR, Multi-Scale Image Registration
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% 
%   - data                 Hands, Omega=(0,20)x(0,25), level=7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - regularizer          mbCurvature
%   - optimization         Gauss-Newton
% ===============================================================================

% Example for Multi-Scale Non-Parametric Image Registration with pre-registration
% (c) Jan Modersitzki 2008/12/30, see FAIRcopyright.m.

% setup data and image viewer, distance, regularizer, pre-registration
setupHandData;  FAIRdiary;
level = 7; omega = MLdata{level}.omega; m = MLdata{level}.m; 
distance('set','distance','SSD');       
regularizer('reset','regularizer','mbElastic','alpha',1e4,'mu',1,'lambda',0);
trafo('reset','trafo','affine2D'); w0 = trafo('w0');

% the y's are used for:  y0/initial guess, yRef/regularization, yStop/stopping
y0   = getStaggeredGrid(omega,m); yRef = y0; yStop = y0;

% discretization of scale space
theta = [logspace(3,0,4),0]; 

% initialize the interpolation scheme and coefficients
inter('set','inter','splineInter','regularizer','moments'); 
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega,'theta',theta(1));
xc    = getCellCenteredGrid(omega,m); 
Rc    = inter(R,omega,xc);

% -- the PIR pre-registration -------------------------
beta = 0; M = []; wRef = []; xc = getCellCenteredGrid(omega,m);
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc); 
[wc,his] = GaussNewton(fctn,w0,'maxIter',500);
reduction = fctn(wc)/fctn(w0);
yc    = grid2grid(trafo(wc,getNodalGrid(omega,m)),m,'nodal','staggered'); 
Yc = {yc}; ITER = max(his.his(:,1)); REDUCTION = reduction;
% parameter for NPIR
NPIRpara = {'maxIter',500,'Plots',@FAIRplots,'yStop',yStop}

% loop over scales
for j=1:length(theta),
  % compute representation of data on j'th scale
  [T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega,'theta',theta(j));
  xc    = getCellCenteredGrid(omega,m); 
  Rc    = inter(R,omega,xc);

  % build objective function and regularizer
  yRef = yc;
  fctn = @(yc) NPIRobjFctn(T,Rc,omega,m,yRef,yc);
  
  % -- solve the optimization problem -------------------------------------------
  FAIRplots('set','mode','NPIR-GN-elastic','omega',omega,'m',m,'fig',j+3,'plots',1);
  FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m));
  [yc,his]  = GaussNewton(fctn,yc,NPIRpara{:});  
  reduction = fctn(yc)/fctn(yStop);
  Yc{end+1} = yc; ITER(end+1) = max(his.his(:,1)); REDUCTION(end+1) = reduction;
end;
diary off