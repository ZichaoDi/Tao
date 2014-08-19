
clear, close all;
viewImage('reset','viewImage','viewImage2D','colormap',gray(256));%,'axis','off');
% ===============================================================================
% dataT1=hdfread('zeolite_4A_2_10x_20mm_00001.hdf','/entry1/data/data','Index',{[1 1],[1 1],[2048 2048]});
% dataT2=hdfread('zeolite_4A_2_10x_20mm_00002.hdf','/entry1/data/data','Index',{[1 1],[1 1],[2048 2048]});
% dataR1=hdfread('zeolite_4A_2_5x_20mm_00001.hdf','/entry1/data/data','Index',{[1 1],[1 1],[2048 2048]});
% dataR2=hdfread('zeolite_4A_2_5x_20mm_00002.hdf','/entry1/data/data','Index',{[1 1],[1 1],[2048 2048]});
% T=double(dataT2)./double(dataT1);
% R=double(dataR2)./double(dataR1);
% ===============================================================================
% T1=double(imread('zeolite_4A_2_5x_00091.tif'));
% T2=double(imread('zeolite_4A_2_5x_00092.tif'));
% T=T2./T1;
% R1=double(imread('zeolite_4A_2_5x_00177.tif'));
% R2=double(imread('zeolite_4A_2_5x_00178.tif'));
% %imagesc(R1),figure,imagesc(R2)
% R=R1./R2;
% return;
% ===============================================================================
global testpad
load CoarsedT256
load CoarsedR256
T=1./CoarsedT256;
R=1./CoarsedR256;
inter('reset','inter','splineInter','regularizer','moments','theta',1e-1);
omega=[0 size(R,1) 0 size(R,2)];
m = floor(size(R)/1);
T = inter('coefficients',T,[],omega,'regularizer','moments','theta',1e-1);
R = inter('coefficients',R,[],omega,'regularizer','moments','theta',1e-1);
T = -log(T);
R=-log(R);
T=GrayScale(T);
R=GrayScale(R);
%  R=T;

testpad=0;%mean(mean(R));
% imagesc(R)
% figure,imagesc(T)
% figure,imshow(mat2gray(T)),figure,imshow(mat2gray(R))
% return;
% ===============================================================================
distance('reset','distance','SSD');
center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rigid2D','c',center);
w0 = trafo('w0');
%  w0= [0 0 100]';
beta = 0;
M =[]; wRef = []; % disable regularization
%M = 2e2*speye(length(w0)); wRef = w0;     % initilize the regularization

% initialize plots
FAIRplots('reset','mode','PIR-GN','fig',1,'scale',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m));
pause;
xc = getCellCenteredGrid(omega,m);
Rc = inter(R,omega,xc);
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc);
% checkDerivative(fctn,w0)
% foo(fctn,w0)
% return;

[wc,his] = GaussNewton(fctn,w0,'Plots',@FAIRplots,'solver','direct','yStop',w0);
% [wc,his] = lBFGS(fctn,w0,'Plots',@FAIRplots,'solver','direct');

