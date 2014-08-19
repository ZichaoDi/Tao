% ===============================================================================
% Example for affine registration 
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
% 
%   - data                 Hands, Omega=(0,20)x(0,25), level=5, m=[32,32]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - optimization         Gauss-Newton
% ===============================================================================

% setup data and initialize image viewer
setupHandData; FAIRdiary; level = 5; omega = MLdata{level}.omega; m = MLdata{level}.m; 

% initialize the interpolation scheme and coefficients
inter('reset','inter','splineInter'); 
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega,'out',0);
xc    = getCellCenteredGrid(omega,m); 
Rc    = inter(R,omega,xc);

% initialize distance measure
distance('set','distance','SSD');       


% -- the PIR pre-registration -------------------------
trafo('reset','trafo','affine2D'); w0 = trafo('w0')
% setup plots and initialize objective function for PIR
FAIRplots('set','mode','PIR-GN-rigid','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 
beta = 0; M = []; wRef = []; xc = getCellCenteredGrid(omega,m);
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc); fctn([]); % report status
% solve the PIR
[wc,his] = GaussNewton(fctn,w0,'maxIter',500,'Plots',@FAIRplots);
