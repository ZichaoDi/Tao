%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: compact version of matrix-free regularization
%
%  S(y) = alpha/2 * norm(B*(y-yRef))^2,
%
% where
%  alpha regularization parameter, weights regularization versus 
%        distance in the joint objective function, alpha = 1 here
%  yRef  is a reference configuration, e.g. yRef = x or a 
%        pre-registration result
%  B     a discrete partial differential operator either in explicit
%        matrix form or as a structure containing the necessary
%        parameters to compute B*y
% see also regularizer E8_regularization_MB
%==============================================================================

clear, close all, help(mfilename);


% initialize the regularization and create a starting  point
regularizer('reset','regularizer','mfElastic','alpha',1,'mu',1,'lambda',0);
y0 = @(omega,m) randn(size(getStaggeredGrid(omega,m)));

% 2D example, initialize physical domain and number of discretization points
omega = [0,1,0,1]; m = [16,12];   % 

% test derivative of 2D implementation
fctn = @(yc) regularizer(yc,omega,m);  checkDerivative(fctn,y0(omega,m));

% 3D example, initialize physical domain and number of discretization points
omega = [0,1,0,1,0,1]; m  = [16,12,8];

% test derivative of 3D implementation
fctn = @(yc) regularizer(yc,omega,m); checkDerivative(fctn,y0(omega,m));