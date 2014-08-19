% ===============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: MLPIR, MultiLevel Parametric Image Registration
% 
%   - data                 Phantom, Omega=(0,200)x(0,200), level=3:7, m=[200,200] 
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     affine2D, not regularized
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all,
viewImage('reset','viewImage','viewImage2D','colormap',gray(256));%,'axis','off'

global TbRs testpad numLevel DrawError
DrawError=0;
numLevel=2;
randn('state',0);

% ===============================================================================

T=double(imread('LoganReference.tiff'));
R=flipud(double(imread('LoganTest1.tiff')));
%  R=double(imread('CroppedRotate.tiff'));
TbRs=1;
% Crop=[100 100 50 50];
% R=double(imcrop(T,Crop));
% omega=[Crop(1)-1 Crop(1)+Crop(3) Crop(2)-1 Crop(2)+Crop(4)];
% omega=20+omega;
inter('reset','inter','splineInter','regularizer','moments','theta',1e0);
Inix=20;
Iniy=30;
omega=[Inix Inix+size(R,1) Iniy Iniy+size(R,2)];
omegat=[0 0+size(T,1) 0 0+size(T,2)];%omega;%
mt=floor(size(T)/1);
m = floor(size(R)/1);
% set view options and interpolation options
viewOptn = {'viewImage','viewImage2D','colormap','bone(256)'};
viewImage('reset',viewOptn{:});
testpad=0;%mean(mean(R));
[MLdata,minLevel,maxLevel,fig] = getMultilevel_phantom({T,R},{omegat,omega},{mt,m});

% % ===============================================================================

distance('reset','distance','SSD');
trafo('reset','trafo','rigid2D');
w0 = trafo('w0');
% run Multilevel Parametric Image Registration
beta = 1; M =[]; wRef=w0; 
% beta = 1e0; M = 5e2*speye(length(w0)); wRef =[0.4,0,0]';     % initilize the regularization
% hd=prod(1./m);
% n=length(w0);
% beta = 1e0; M = 5e3*spdiags(ones(n,1)*[-1,2,-1],-1:1,n,n); wRef = [0.3,0,0]'; 
wc = MLPIR_phantom(MLdata,minLevel,maxLevel,'M',M,'beta',beta,'wRef',wRef,'plotIter',0,'plotMLiter',1);



