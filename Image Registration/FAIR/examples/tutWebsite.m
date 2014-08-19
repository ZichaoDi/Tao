% Tutorial for FAIR: All in One
% (c) Jan Modersitzki 2009/03/24, see FAIR.2 and FAIRcopyright.m.
% - setup data [ShekharEtAl2005]
% - generate multi-level representation  (see getMultilevel.m)
% - setup  viewer          (viewImage2D), 
%           interpolator    (splineInter, regularized), 
%           distance        (NGF), 
%           transformation  (affine2D, not regularized)
%           regularizer     (curvature, matrix free)
% - run Multi-Level Iage Registration 
clc, clear; close all; help(mfilename)

%% -- start getting the data
% reads and converts data
image = @(str) double(flipud(imread(str))'); % 

% Original data from [ShekharEtAl2005]
%@article{ShekharEtAl2005,
%   author = {Raj Shekhar and Vivek Walimbe and Shanker Raja and Vladimir Zagrodsky
%             and Mangesh Kanvinde and Guiyun Wu and Bohdan Bybel},
%    title = {Automated 3-Dimensional Elastic Registration of Whole-Body {PET}
%             and {CT} from Separate or Combined Scanners},
%  journal = {J. of Nuclear Medicine},
%   volume = {46},
%   number = {9},
%     year = {2005},
%    pages = {1488--1496},
%}
% load the original data, convert from uint8 -> double, ij -> xy
dataT = double(flipud(imread('PET-CT-PET.jpg'))');
dataR = double(flipud(imread('PET-CT-CT.jpg'))');

%set physical domain, initialize discretization and grid
omega = [0,50,0,50];
m     = [128,128];
xc    = getCellCenteredGrid(omega,m);

% setup image viewer  
viewImage('reset','viewImage','viewImage2D','colormap','bone(256)');

% visualize data
figure(1); clf;
subplot(1,2,1); viewImage(dataT,omega,size(dataT));
th(1) = gca; th(2) = title('T=PET');
subplot(1,2,2); viewImage(dataR,omega,size(dataR)); hold on;
th(3) = gca; th(4) = title('R=CT'); 
plotGrid(xc,omega,m,'spacing',[8,8]);
set(th,'fontsize',30);

% create multilevel representation of the data
MLdata = getMultilevel({dataT,dataR},omega,m);

%  -- end getting the data

%% -- initialize interpolation, distance, transformation, regularization
%     lots to play with!
theta = 5e-3; % play with it
inter('reset','inter','splineInter','regularizer','moments','theta',theta);
distance('reset','distance','NGF'); 
trafo('reset','trafo','affine2D');    
regularizer('reset','regularizer','mfElastic','alpha',1e-1,'mu',1','lambda',0');
% --  end initialization

%% -- Multi-Level Image Registration
yc =  MLIR(MLdata,'minLevel',4,'plotIter',0,'plotMLiter',0);

%% -- visualize results
showResults(MLdata,yc,'fig',10);

