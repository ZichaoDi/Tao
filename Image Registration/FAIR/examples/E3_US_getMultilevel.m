%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: creates a multilevel data representation
%
%==============================================================================

clear, close all, help(mfilename); 
echo on

% load data, define a doman and an initial discretization
dataT  = double(imread('US.jpg'));
omega  = [0,size(dataT,1),0,size(dataT,2)];
m      = [128,128];

% setup image viewer
viewImage('reset','viewImage','viewImage2D','colormap','gray(256)');

% creating a multilevel representation using getMultilevel.m
MLdata = getMultilevel(dataT,omega,m);

% display multilevel representation of level 3
disp('MLdata{3}='); disp(MLdata{3})
echo off