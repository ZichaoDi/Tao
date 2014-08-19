%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: curvature regularization
%
% illustrates the tools for L2-norm based regularization
%
%  S(y) = alpha/2 * norm(B*(y-yRef)^2,
%
% where
%  alpha regularization parameter, weights regularization versus 
%        distance in the joint objective function, alpha = 1 here
%  yRef  is a reference configuration, e.g. yRef = x or a 
%        pre-registration result
%  B     a discrete partial differential operator either in explicit
%        matrix form or as a structure containing the necessary
%        parameters to compute B*y
% see also regularizer and getCurvatureMatrix
%==============================================================================

clear, close all, help(mfilename);

alpha  = 1; % regularization parameter, irrelevant in this tutorial

% 3D example
omega  = [0,4,0,2,0,1]; % physical domin
m      = [16,12,8];     % number of discretization points
yRef   = getCellCenteredGrid(omega,m); % reference for regularization
yc     = rand(size(yRef));          % random transformation

% build elasticity operator on a staggered grid
B = getCurvatureMatrix(omega,m);
FAIRfigure(1); clf; subplot(1,3,1); spy(B); title('B curvature on cell-centered grid')


regularizer('reset','regularizer','mbCurvature','alpha',alpha);
regularizer('disp');

hd  = prod((omega(2:2:end)-omega(1:2:end))./m);
mbA = hd*alpha*(B'*B);

[S,dS,d2S] = regularizer(yc-yRef,omega,m);
subplot(1,3,1); spy(mbA); title('hd\alpha B''*B')
subplot(1,3,2); spy(d2S); title('d2S')
subplot(1,3,3); spy(mbA-d2S); title('difference')

Splain = 0.5*(yc-yRef)'*mbA*(yc-yRef)
Smb    = regularizer(yc-yRef,omega,m)
regularizer('set','regularizer','mbCurvature');
Smf    = regularizer(yc-yRef,omega,m)

