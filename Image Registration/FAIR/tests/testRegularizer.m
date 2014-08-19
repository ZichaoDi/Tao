% (c) Jan Modersitzki 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Test regualarization toolbox

toolbox = fileparts(which('regularizer'));
testStart

debit = {
          'contents.m'...
          'curvature.m'...
          'elastic.m'...
          'geometry.m'...
          'geometrymexC.cpp'...
          'geometrymexC.h'...
          'getCurvatureMatrix.m'...
          'getElasticMatrixNodal.m'...
          'getElasticMatrixStg.m'...
          'getGradientNodal.m'...
          'hyperElastic.m'...
          'regularizer.m'...
  };

checkToolbox(toolbox,debit,'EDIT',0,'KEYBOARD',0,'RUN',0,'MEX',0)

%% parameters
% used in all regularizers
alpha  = 1;

% for (linear) elastic regularizer
mu     = 1;
lambda = 0;

% for hyperelastic regularizer
alphaLength = 13.2;
alphaArea   = 1.2;
alphaVolume = 3.2;


%% test syntax of regularizer
% linear elastic
regularizer('clear');
regularizer('reset','regularizer','mfElastic',...
  'alpha',alpha,'mu',mu,'lambda',lambda);
regularizer('disp');
[scheme,parameter] = regularizer;
parameter = cell2struct({parameter{2:2:end}},{parameter{1:2:end}},2)

% curvature
regularizer('clear');
regularizer('reset','regularizer','mbCurvature','alpha',alpha);
regularizer('disp');
[scheme,parameter] = regularizer;
parameter = cell2struct({parameter{2:2:end}},{parameter{1:2:end}},2)

% hyperelastic
regularizer('clear');
regularizer('reset','regularizer','mfHyperElastic','alpha',alpha,...
    'alphaLength',alphaLength,'alphaArea',alphaArea,'alphaVolume',alphaVolume);
regularizer('disp');
[scheme,parameter] = regularizer;
parameter = cell2struct({parameter{2:2:end}},{parameter{1:2:end}},2)

%%
omega = [0,1,0,1];
m = [32,32]/8;
hd = prod((omega(2:2:end)-omega(1:2:end))./m);
regularizer('reset','regularizer','mfElastic',...
  'alpha',alpha,'mu',mu,'lambda',lambda);
% check elasticity
H.omega  = omega; 
H.m      = m; 
H.alpha  = regularizer('get','alpha');
H.mu     = regularizer('get','mu');
H.lambda = regularizer('get','lambda');
H.regularizer = regularizer;


X   = getStaggeredGrid(omega,m);
B  = getElasticMatrixStg(omega,m,mu,lambda);
Y  = randn(size(X));

M   = spdiags(1+0*Y,0,length(Y),length(Y));
A   = M + hd* H.alpha *B'*B;
rhs = A*X;
% rhs = A(:,5);
u0  = zeros(size(Y));

% prepare for multigrid
H.MGlevel      = log2(m(1))+1;
H.MGcycle      = 4;
H.MGomega      = 2/3;   %% !!! 0.5 should be better
H.MGsmoother   = 'mfJacobi';
H.MGpresmooth  = 10;
H.MGpostsmooth = 10;
H.d2D.M        = full(diag(M));

[Sc, dS, d2S] = regularizer(X,omega,m);
H.d2S = d2S;

Zmg = mfvcycle(H,u0,rhs,1e-12,H.MGlevel,5);
testMG = norm( A\rhs-Zmg);
fprintf('|A\\rhs-uMG|=%s\n',num2str(testMG));


%% check matrix free implementation
Omega = {[0,1,0,2],[0,1,0,2,0,3]};
M = {[4,5],[3,4,5]};
alpha = 1; mu = 1; lambda = 0;

for j=1:length(Omega),
  
  omega = Omega{j}; m = M{j}; 
  fprintf('\ndimension d=%d\n',length(omega)/2)
  
  
  % elastic
  fprintf('the elastic version\n')
  X = getStaggeredGrid(omega,m);
  Y = randn(size(X));
  regularizer('reset','regularizer','mfElastic','alpha',alpha,'mu',mu,'lambda',lambda);
  [mfS,mfdS,mfd2S] = regularizer(Y,Omega{j},M{j});
  
  regularizer('reset','regularizer','mbElastic','alpha',alpha,'mu',mu,'lambda',lambda);
  [mbS,mbdS,mbd2S] = regularizer(Y,Omega{j},M{j});
 
  testS = norm(mfS - mbS) / norm(mbS);
  fprintf('%-30s = %s\n','norm(mfS - mbS) / norm(mbS)',num2str(testS))
  
  testdS = norm(mfdS - mbdS) / norm(mbdS);
  fprintf('%-30s = %s\n','norm(mfdS - mbdS) / norm(mbdS)',num2str(testdS))
  
  % curvature 
  fprintf('the curvature version\n')
  X = getCellCenteredGrid(omega,m);
  Y = randn(size(X));
  regularizer('reset','regularizer','mfCurvature','alpha',alpha);
  [mfS,mfdS,mfd2S] = regularizer(Y,Omega{j},M{j});

  
  regularizer('reset','regularizer','mbCurvature','alpha',alpha);
  [mbS,mbdS,mbd2S] = regularizer(Y,Omega{j},M{j});
   
  testS = norm(mfS - mbS) / norm(mbS);
  fprintf('%-30s = %s\n','norm(mfS - mbS) / norm(mbS)',num2str(testS))
  
  testdS = norm(mfdS - mbdS) / norm(mbdS);
  fprintf('%-30s = %s\n','norm(mfdS - mbdS) / norm(mbdS)',num2str(testdS))

  % hyperelastic
  fprintf('the hyperelastic version\n')
  X = getNodalGrid(omega,m);
  Y = randn(size(X));
  regularizer('reset','regularizer','mfHyperElastic','alpha',alpha,...
    'alphaLength',alphaLength,'alphaArea',alphaArea,'alphaVolume',alphaVolume);
  [mfS,mfdS,mfd2S] = regularizer(Y,Omega{j},M{j});

  
   regularizer('reset','regularizer','mbHyperElastic','alpha',alpha,...
    'alphaLength',alphaLength,'alphaArea',alphaArea,'alphaVolume',alphaVolume);
   [mbS,mbdS,mbd2S] = regularizer(Y,Omega{j},M{j});
   
  testS = norm(mfS - mbS) / norm(mbS);
  fprintf('%-30s = %s\n','norm(mfS - mbS) / norm(mbS)',num2str(testS))
  
  testdS = norm(mfdS - mbdS) / norm(mbdS);
  fprintf('%-30s = %s\n','norm(mfdS - mbdS) / norm(mbdS)',num2str(testdS))

end;
%% check derivatives
Omega = {[0,1,0,2],[0,1,0,2,0,3]};
M = {[4,5],[3,4,5]};
alpha = 1; mu = 1; lambda = 0;

for j=1:length(Omega),
  
  omega = Omega{j}; m = M{j}; 
  dim   = length(omega)/2;
  fprintf('\ndimension d=%d\n',dim)
  
  % elastic
  fprintf('the elastic version\n')
  X = getStaggeredGrid(omega,m);
  Y = randn(size(X));
  regularizer('reset','regularizer','mfElastic','alpha',alpha,'mu',mu,'lambda',lambda);
  fctn = @(Y) regularizer(Y,Omega{j},M{j});
  checkDerivative(fctn,Y);
  title(sprintf('checkDerivative: elastic - %dD',dim));
  
  % curvature
  fprintf('the curvature version\n')
  X = getCellCenteredGrid(omega,m);
  Y = randn(size(X));
  regularizer('reset','regularizer','mfCurvature','alpha',alpha);
  fctn = @(Y) regularizer(Y,Omega{j},M{j});
  checkDerivative(fctn,Y);
  title(sprintf('checkDerivative: curvature - %dD',dim));
  
  % hyperelastic
  fprintf('the hyperelastic version\n')
  X = getNodalGrid(omega,m);
  Y = randn(size(X));
  regularizer('reset','regularizer','mfHyperElastic','alpha',alpha,...
    'alphaLength',alphaLength,'alphaArea',alphaArea,'alphaVolume',alphaVolume);
  fctn = @(Y) regularizer(Y,Omega{j},M{j});
  checkDerivative(fctn,Y);
  title(sprintf('checkDerivative: hyperelastic - %dD',dim));
  
end;

testEnd

%{ 
	=======================================================================================
	FAIR: Flexible Algorithms for Image Registration, Version 2011
	Copyright (c): Jan Modersitzki
	Maria-Goeppert-Str. 1a, D-23562 Luebeck, Germany
	Email: jan.modersitzki@mic.uni-luebeck.de
	URL:   http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
	=======================================================================================
	No part of this code may be reproduced, stored in a retrieval system,
	translated, transcribed, transmitted, or distributed in any form
	or by any means, means, manual, electric, electronic, electro-magnetic,
	mechanical, chemical, optical, photocopying, recording, or otherwise,
	without the prior explicit written permission of the authors or their
	designated proxies. In no event shall the above copyright notice be
	removed or altered in any way.

	This code is provided "as is", without any warranty of any kind, either
	expressed or implied, including but not limited to, any implied warranty
	of merchantibility or fitness for any purpose. In no event will any party
	who distributed the code be liable for damages or for any claim(s) by
	any other party, including but not limited to, any lost profits, lost
	monies, lost data or data rendered inaccurate, losses sustained by
	third parties, or any other special, incidental or consequential damages
	arrising out of the use or inability to use the program, even if the
	possibility of such damages has been advised against. The entire risk
	as to the quality, the performace, and the fitness of the program for any
	particular purpose lies with the party using the code.
	=======================================================================================
	Any use of this code constitutes acceptance of the terms of the above statements
	=======================================================================================
%}