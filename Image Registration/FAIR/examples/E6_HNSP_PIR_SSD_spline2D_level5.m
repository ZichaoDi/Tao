% ===============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: Parametric Image Registration
% 
%   - data                 HNSP, Omega=(0,2)x(0,1), level=5, m=[256,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     splineTransformation2D
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);

setupHNSPData; 
inter('set','inter','splineInter'); 
level = 5; omega = MLdata{level}.omega; m = MLdata{level}.m; 
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
xc = getCellCenteredGrid(omega,m); 
Rc = inter(R,omega,xc);

distance('set','distance','SSD');       % initialize distance measure

% initialize transformation and starting guess; 
% here: spline with 2 times [4,5] coefficients
trafo('reset','trafo','splineTransformation2D','omega',omega,'m',m,'p',[4,5]);
w0 = trafo('w0'); 

FAIRplots('reset','mode','PIR-spline','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% ----- call Gaus-Newton ------------------------------------
GNoptn = {'maxIter',50,'Plots',@FAIRplots};
fctn  = @(wc) PIRobjFctn(T,Rc,omega,m,0,[],[],xc,wc);
[wc,his] = GaussNewton(fctn,w0,GNoptn{:});

figure(1); clf
viewImage(inter(T,omega,xc),omega,m,'axis','off'); hold on
ph = plotGrid(trafo(wc,xc),omega,m,'spacing',1,'linewidth',1,'color','w');
