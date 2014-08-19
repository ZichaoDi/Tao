%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: demo for 2D data setup and visualization
%
% load images                        (hands-?.jpg)
% visualize ij and omega based data  (viewImage)
% creates multi-lebel representation (getMultilevel) 
%==============================================================================

clear, close all, help(mfilename); echo on

% load data
Tij = imread('hands-T.jpg'); Tdata = double(flipud(Tij))';
Rij = imread('hands-R.jpg'); Rdata = double(flipud(Rij))';

% specify domain Omega = [omega(1),omega(2)]x[omega(3),omega(3)];
omega = [0,20,0,25];
m     = size(Tdata);

% setup image viewer
viewImage('reset','viewImage','viewImage2D','colormap','gray(256)');

% visualize
FAIRfigure(1); colormap(gray(256));
subplot(2,2,1); imagesc(Tij);              title('original T data, uint8, ij');
subplot(2,2,2); viewImage(Tdata,omega,m);  title('FAIR T data on \Omega, xy');
subplot(2,2,3); imagesc(Rij);              title('original R data, uint8, ij');
subplot(2,2,4); viewImage(Rdata,omega,m);  title('FAIR R data on \Omega, xy');

% create multi-level representaion
MLdata = getMultilevel({Tdata,Rdata},omega,m,'fig',2);

echo off
