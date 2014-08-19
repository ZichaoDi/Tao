% ===============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: MLPIR, Parametric Image Registration, scale space
% 
%   - data                 HNSP, Omega=(0,2)x(0,1), level=5, m=[256,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - optimization         GaussNewton
% ===============================================================================

clear, close all, help(mfilename);

setupHNSPData; 
inter('set','inter','splineInter'); 
level = 5; omega = MLdata{level}.omega; m = MLdata{level}.m; 
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
xc = getCellCenteredGrid(omega,m); 
Rc = inter(R,omega,xc);

distance('set','distance','SSD');       % initialize distance measure

trafo('reset','trafo','affine2D');
w0 = trafo('w0'); 

FAIRplots('set','mode','PIR-affine','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% -- solve the optimization problem on different scales -----------------------

wStop = w0;             % use one GLOBAL stopping criterion
theta = [1e2,1e1,1,0];  % the different scales: smooth to detailed
for j=1:length(theta),
  % initilaize the data for scale theta(j)
  inter('reset','inter','splineInter','regularizer','moments','theta',theta(j));
  [T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
  
  % setup plots and initialize it
  FAIRplots('reset','mode','PIR-scale','fig',j,'plots',1);
  FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m));
  
  % initialize optimizer
  Rc   = inter(R,omega,xc); % note Rc depends on the scale as well
  fctn = @(wc) PIRobjFctn(T,Rc,omega,m,0,[],[],xc,wc); 
  if j ==1, fctn([]); end;  % report status
  
  % solve problem for this scale 
  [wc,his] = GaussNewton(fctn,w0,'yStop',wStop,'Plots',@FAIRplots);
  % and use solutions as starting guess for next scale  
  w0 = wc;
end;
