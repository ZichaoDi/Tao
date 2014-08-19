% (c) Jan Modersitzki 2010/12/25, see FAIR.2 and FAIRcopyright.m.
% Flexible Algorithms for Image Registration, SIAM, 2009
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% Loads ultrasound image.
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

example = 'US';
checkSetupDataFile; if OK, return; end;

txt = {
  'creates '
  ' - dataT:    data for templateimage, d-array of size p'
  ' - omega:    coding domain = [omega(1),omega(2)]x...x[omega(2*d-1),omega(2*d)]'
  ' - m:        size of the finest interpolation grid'
  ' - MLdata:   MLdata representation of the data'
  };
fprintf('%s\n',txt{:});

% do whatever needed to be done to get your data here
image = @(str) double(flipud(imread(str))); % reads and converts

% load the original data, set domain, initial discretization, and grid
dataT = image('US.jpg');
omega = [0,size(dataT,1),0,size(dataT,2)];
m     = 128*[3,2];

% set view options and interpolation options
viewOptn = {'viewImage','viewImage2D','colormap','bone(256)'};
viewImage('reset',viewOptn{:});
intOptn = {'inter','linearInter'};
inter('reset',intOptn{:});

str = 'Vibe Heldmann.jpg';
viewData  = @(I) viewImage(inter(I,omega,getCellCenteredGrid(omega,m)),omega,m);

FAIRfigure(1,'figname',mfilename); clf;
viewData(dataT); title(str,'interpreter','none');
MLdata = getMultilevel(dataT,omega,m,'fig',2);

% save to outfile
save(outfile,'dataT','omega','m','MLdata','viewOptn','intOptn');
