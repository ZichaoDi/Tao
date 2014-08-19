% (c) Jan Modersitzki 2010/12/25, see FAIR.2 and FAIRcopyright.m.
% Flexible Algorithms for Image Registration, SIAM, 2009
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% This extended example indicates the data format used in FAIR.
%
% This file adds the following variables to outfile.mat:
%   dataT      data of template  image
%   dataR      data of reference image
%   omega      = [left,right,bottom,top]  domain description
%   m          size for finest data representation
%   MLdata     multilevel representation of data, see getMultilevel for details
%   viewOptn   parameterization of image viewer
%   intOptn  parameterization of image interpolater
%
% see also get2Ddata and data/contents.m

txt = {
  'creates '
  ' - dataT:    data for template image, d-array of size p'
  ' - dataR:    data for reference image, d-array of size q'
  ' - omega:    coding domain = [omega(1),omega(2)]x...x[omega(2*d-1),omega(2*d)]'
  ' - m:        size of the finest interpolation grid'
  ' - MLdata:   MLdata representation of the data'
  };
fprintf('%s\n',txt{:});

outfile = fullfile(FAIRpath,'temp',[mfilename,'.mat']);
example = 'PET-CT';

% here, original data from
% 	@article{ShekharEtAl2005,
% 	 author = {Raj Shekhar and Vivek Walimbe and Shanker Raja and Vladimir Zagrodsky
% 	           and Mangesh Kanvinde and Guiyun Wu and Bohdan Bybel},
% 	 title = {Automated 3-Dimensional Elastic Registration of Whole-Body {PET}
% 	          and {CT} from Separate or Combined Scanners},
% 	 journal = {J. of Nuclear Medicine},
% 	 volume = {46},
% 	 number = {9},
% 	 year = {2005},
% 	 pages = {1488--1496},
% 	}

% do whatever needed to be done to get your data here
image = @(str) double(flipud(imread(str))'); % reads and converts

% load the original data, set domain, initial discretization, and grid
dataT = image('PET-CT-PET.jpg');
dataR = image('PET-CT-CT.jpg');
omega = [0,50,0,50]; % [left,right,bottom,top]
m     = [128,128];

% set view options and interpolation options
viewOptn = {'viewImage','viewImage2D','colormap','bone(256)'};
viewImage('reset',viewOptn{:});

intOptn = {'inter','linearInter'};
inter('reset',intOptn{:});

% some plot
xc = getCellCenteredGrid(omega,m);
viewData  = @(I,omega) viewImage(inter(I,omega,xc),omega,m);

FAIRfigure(1,'figname',mfilename); clf;
subplot(1,2,1); viewData(dataT,omega); title('template');
subplot(1,2,2); viewData(dataR,omega); title('reference');

% create multilevel representation of the data
MLdata = getMultilevel({dataT,dataR},omega,m,'fig',2);

% save to outfile
save(outfile,'dataT','dataR','omega','m','MLdata','viewOptn','intOptn');