%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: shows how to visulaize 2D data with FAIR.
%
% - load data             ('US.jpg')
% - setup viewer          (viewImage2D)
% - view  data
%==============================================================================

clear, close all, help(mfilename);
echo on

% load data
dataT = double(imread('US.jpg'));
m     = size(dataT);
omega = [0,m(1),0,m(2)]; 

% setup image viewer
viewImage('reset','viewImage','viewImage2D','colormap','gray(256)','axis','off');

% view data
viewImage(dataT,omega,m,'title','my first image');

echo off