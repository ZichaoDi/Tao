%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: Initializes 3D data and show how this can be used in FAIR.
%
% This file creates an the following data (which can then be saved): 
%   dataT     template  image, a d-array of size m, 
%   dataR     reference image, a d-array of size m, 
%   omega     domain description 
%   m         initial discretization 
%   ML        multi-level representation of the data
%
%   viewOptn  options for image viewer
%   intOptn options for image interpolation
%==============================================================================

% load 3D data
load brain3D; 
whos
viewOptn  = {'viewImage','imgmontage'};  viewImage('reset',viewOptn{:});
intOptn = {'inter','linearInter'};   inter('reset',intOptn{:});
figure(2); colormap(gray(256));
MLdata    = getMultilevel({dataT,dataR},omega,m);
save brain3DML.mat dataT dataR omega m MLdata viewOptn intOptn

