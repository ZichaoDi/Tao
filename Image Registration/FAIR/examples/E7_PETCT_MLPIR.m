%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: various distances and Multi-Level Parametric Image Registration
%
%   - data                 PETCT, Omega=(0,140)x(0,151), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             {'SSD','NCC','MI','NGF'}
%   - transformation       affine2D
% see also E7_PETCT_MLPIR_ext
%==============================================================================

clear, close all, help(mfilename);
setupPETCTData;                                         % load data

% a list of distance measures to be used
DM = {'SSD','NCC','MI','NGF'};

for dm = 1:length(DM), % run over all distance measures

  % initialize interpolation, using a smooth representation (theta=1e0) 
  inter('reset','inter','splineInter','regularizer','moments','theta',1e0);

  % initialize transformation, create initial guess and reference for stopping
  trafo('reset','trafo','affine2D'); wStop = trafo('w0'); w0 = wStop;

  % initialize distance and display options
  distance('reset','distance',DM{dm}); distance('disp')
  
  % run MLPIR using sufficient amount of details (level=5)
  wSmooth =  MLPIR(MLdata,'minLevel',5,'plotIter',0,'plotMLiter',0);

  % refine interpolation (theta=1e-3)
  inter('reset','inter','splineInter','regularizer','moments','theta',1e-3);
  level = length(MLdata); omega = MLdata{level}.omega; m = MLdata{level}.m;
  [T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
  
  % start PIR, using the result from the smooth problem as starting guess
  
  % initialize plots
  FAIRplots('set','mode','PIR');
  FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m));

  % optimize
  xc = getCellCenteredGrid(omega,m);   
  Rc = inter(R,omega,xc);
  fctn = @(wc) PIRobjFctn(T,Rc,omega,m,0,[],[],xc,wc); fctn([]);
  [wc,his] = GaussNewton(fctn,w0,'Plots',@FAIRplots);

  % visualize results
  yc = trafo(wc,xc);
  R0 = inter(R,omega,xc);
  T0 = inter(T,omega,xc);
  Tc = inter(T,omega,yc);

  figure(11); clf;
  viewImage(T0,omega,m,'axis','off'); hold on;
  plotGrid(yc,omega,m,'spacing',ceil(m/32),'linewidth',2,'color','w');

  figure(12); clf;
  overlayImage2D(Tc,R0,omega,m); axis off;
end;
