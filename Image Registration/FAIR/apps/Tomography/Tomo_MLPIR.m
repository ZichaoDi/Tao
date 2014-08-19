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

clear, close all, help(mfilename);

% load data, initialize image viewer, interpolator, transformation, distance
T1=double(imread('zeolite_4A_2_5x_00091.tif'));
T2=double(imread('zeolite_4A_2_5x_00092.tif'));
T=T2./T1;
R1=double(imread('zeolite_4A_2_5x_00177.tif'));
R2=double(imread('zeolite_4A_2_5x_00178.tif'));
R=R2./R1;
%%%=========================================================
omega=[0 size(R,1) 0 size(R,2)];
m = floor(size(R)/1);
% set view options and interpolation options
viewOptn = {'viewImage','viewImage2D','colormap','bone(256)'};
viewImage('reset',viewOptn{:});

[MLdata,minLevel,maxLevel,fig] = getMultilevel({T,R},omega,m);
return; 
% % ===============================================================================
% T=double(imread('LenaCropped.tiff'));
% R=double(imread('LenaReference.tiff'));
% omega=[0 size(R,1) 0 size(R,2)];
% m = floor(size(R)/1);
% % ===============================================================================
inter('reset','inter','linearInter','regularizer','moments','theta',1e-1);
distance('reset','distance','SSD');
trafo('reset','trafo','rigid2D');
w0 = trafo('w0');
% run Multilevel Parametric Image Registration
beta = 1; M =[];% 5e2*speye(numel(w0)); 
wRef = w0;
wc = MLPIR(MLdata,'M',M,'beta',beta,'wRef',wRef,'plotIter',0,'plotMLiter',1);