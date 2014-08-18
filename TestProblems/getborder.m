function [x,y,hx,hy,sx0,sy0,sx1,sy1] = getborder(nx,ny)
%--------------------------------------------------------------
% specify boundary conditions for optimal control problem
%
% Parameters:
%
% Input:
%	nx	- # of points in x direction
%	ny	- # of points in y direction
% Output:
%	x	- array of x values
%	y	- array of y values
%	hx	- spacing in x direction
%	hy	- spacing in y direction

x0 = 0;
x1 = 1;
y0 = 0;
y1 = 1;
hx = (x1 - x0)/(nx-1);
hy = (y1 - y0)/(ny-1);
x1  = (x0:hx:x1)';
y1  =(y0:hy:y1);
x=x1(2:end-1);
y=y1(2:end-1);
sx0=0*x1';
sy0=0*y1';
sx1=0*x1';
sy1=0*y1';
%-------------------------------------------------------
