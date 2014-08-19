

clear, close all,
viewImage('reset','viewImage','viewImage2D','colormap',gray(256));%,'axis','off'
% ===============================================================================
global testpad DrawError regEnable
testpad=0;%mean(mean(R));
DrawError=1;
randn('state',0);
c=[3 7]*1;
rec=[10 10]*1;
r=3;
[T,R]=CreateSimpleTests(c,rec,r);
% initialize the interpolation scheme and coefficients 
inter('reset','inter','linearInter','regularizer','moments','theta',1e-1);
omega=[0 size(R,1) 0 size(R,2)];
omegat=omega;%[10 10+size(T,1) 40 40+size(T,2)];%[Crop(1)-1 Crop(1)+Crop(3) Crop(2)-1 Crop(2)+Crop(4)];%
mt=floor(size(T)/1);
m = floor(size(R)/1);
T = inter('coefficients',T,[],omegat,'regularizer','moments','theta',1e-1);
R = inter('coefficients',R,[],omega,'regularizer','moments','theta',1e-1);
xc    = getCellCenteredGrid(omega,m); 
Rc    = inter(R,omega,xc);

% initialize distance measure
distance('set','distance','SSD');       

% initialize regularization, note: yc-yRef is regularized, elastic is staggered 
regEnable=1;
regularizer('reset','regularizer','mbCurvature','alpha',1e4,'mu',1,'lambda',0);
% y0   = getStaggeredGrid(omega-0.2,m);
y0   = getCellCenteredGrid(omega-1,m); 
yRef = y0; yStop = y0;


% setup and initialize plots 
FAIRplots('set','mode','NPIR-Gauss-Newton','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% build objective function, note: T coefficients of template, Rc sampled reference
fctn = @(yc) NPIRobjFctn(T,Rc,omega,m,yRef,yc); fctn([]); % report status

% -- solve the optimization problem -------------------------------------------
[yc,his] = GaussNewton(fctn,y0,'maxIter',500,'Plots',@FAIRplots,'yStop',yStop);
% report results
% iter = size(his.his,1)-2; reduction = 100*fctn(yc)/fctn(y0);
% fprintf('reduction = %s%% after %d iterations\n',num2str(reduction),iter);
