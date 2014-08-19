% ==================================================================================
% (c) Jan Modersitzki 2010/12/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function [Sc,dS,d2S] = curvature(uc,omega,m,varargin)
%
% computes curvature regularization energy for uc = yc - yRef
% (cell centered)
%
% S(u) = 0.5 * \int_{\omega} u(x)' * B' * B * u(x) dx,
%
% where B is the curvature operator, see getCurvatureMatrix .
%
% Input:
% ------
%   uc       	displacement field (staggered)
%   omega       spatial domain
%   m     		number of discretization points
%   varargin	optional parameters (see below)
%
% Output:
% -------
%   Sc          current value  (0.5 * hd * uc'*B'*B*uc)
%   dS          derivative     (hd * uc'*B'*B)
%   d2S   		Hessian        (B'*B)
%  if ~matrixFree,  d2S is sparse matrix; else, d2S is struct endif		
%
% see also BigTutorialRegularizer
% ==================================================================================
function [Sc,dS,d2S] = curvature(uc,omega,m,varargin)
if nargin == 0
    help(mfilename);
    runMinimalExample;
    return;
end
persistent A omegaOld mOld

if ~exist('mOld','var'),     mOld = [];     end;
if ~exist('omegaOld','var'), omegaOld = []; end;

matrixFree  = 0;
alpha       = 1;

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

hd  = prod((omega(2:2:end)-omega(1:2:end))./m);

if not(matrixFree), % matrix-based
  build = isempty(mOld) || isempty(omegaOld) ...
    || length(mOld) ~= length(m) || length(omegaOld) ~= length(omega) ...
    || any(mOld ~= m) || any(omegaOld~=omega);
  if build,
    mOld = m; omegaOld = omega;
    A = getCurvatureMatrix(omega,m);
    A = alpha*hd*(A'*A);
  end;
  dS  = uc'*A;
  Sc  = 0.5*dS*uc;
  d2S = A;
else % matrix-free
  d2S.regularizer = regularizer;
  d2S.alpha  = alpha;
  d2S.By     = @(u,omega,m) curvatureOperator(u,omega,m,'By');
  d2S.BTy    = @(u,omega,m) curvatureOperator(u,omega,m,'BTy');
  d2S.B      = @(omega,m)getCurvatureMatrix(omega,m);
  d2S.diag   = @(omega,m) getDiag(omega,m,alpha);
  d2S.solver = 'PCG-curvature';

  d2S.d2S  = @(uc,omega,m) ...
             alpha * prod((omega(2:2:end)-omega(1:2:end))./m) * ...
             curvatureOperator(curvatureOperator(uc,omega,m,'By'),omega,m,'BTy');
 
  
  dS   = d2S.d2S(uc,omega,m)';
  Sc   = .5*dS*uc;
end

function D = getDiag(omega,m,alpha)
h    = (omega(2:2:end)-omega(1:2:end))./m;
hd   = prod(h);
dim  = length(omega)/2;
D    = (2/h(1)^2+2/h(2)^2)^2*ones(size(prod(m)*dim));
D    = alpha*hd*spdiags(D,0,numel(D),numel(D));

function By = curvatureOperator(uc,omega,m,flag)
dim = length(omega)/2;
flag = sprintf('%s-%dD',flag,dim);
switch flag,
  case {'By-2D','BTy-2D'},
    %
    % By = | \Delta y^1 0          | = | d_1^2 y^1 + d_2^2 y^1 |
    %      | 0          \Delta y^2 | = | d_1^2 y^2 + d_2^2 y^2 |   
    uc = reshape(uc,[m,2]);
    By = zeros(size(uc));
    d2 = @(i) D2(i,omega,m);
    By(:,:,1) = d2(1)*uc(:,:,1) + uc(:,:,1)*d2(2);
    By(:,:,2) = d2(1)*uc(:,:,2) + uc(:,:,2)*d2(2);
    By = By(:);
  
  case {'By-3D','BTy-3D'},
    %
    %      | \Delta y^1 0          0          | = | d_1^2 y^1 + d_2^2 y^1 + d_3^2 y^1 |
    % By = | 0          \Delta y^2 0          | = | d_1^2 y^2 + d_2^2 y^2 + d_3^2 y^2 |   
    %      | 0          0          \Delta y^3 | = | d_1^2 y^3 + d_2^2 y^3 + d_3^2 y^3 |   

    % efficient implementation of Bcurvature*Y	
    n  = prod(m);            % number of voxels
    uc = reshape(uc,[m,3]);  % note that now uc(:,:,:,ell) is the ell-th
                             % component of Y=(Y^1,Y^2,Y^3)
                             % of size m(1)-by-m(2)-by-m(3)
    By = zeros(numel(uc),1); % allocate memory for the output
    d2 = @(i) D2(i,omega,m); % this is a shortcut to the discrete second derivative

    % the following line is a shortcut for 
    %  - permuting the 3d-array using the permutation J
    %  - reshape it to a 2D-array of size q-by-prod(m)/q, where q=m(J(1))
    %  - multiply by A (which is q-by-q)
    %  - undo the reshape, i.e. make the result to m(J(1))-by-m(J(2))-by-m(J(3))
    %  - undo the permutation
    operate = @(A,z,J) ipermute(reshape(A*reshape(permute(z,J),m(J(1)),[]),m(J)),J);

    % run over all compunents y^ell of Y=(Y^1,Y^2,Y^3)
    for ell=1:3,
      % compute
      % (I_3\otimes I_2\otimes d2(1) + I_3\otimes d2(2)\otimes I_1 ...
      % + d2(3)\otimes I_2\otimes I_1) y^ell
      for k=1:3,
        z = operate(d2(k),uc(:,:,:,ell),[k,setdiff(1:3,k)]);
        By((ell-1)*n+(1:n)) = By((ell-1)*n+(1:n)) + reshape(z,[],1);
      end;
    end;   
end;

function D = D2(i,omega,m)
h = (omega(2:2:end)-omega(1:2:end))./m;
D = spdiags(ones(m(i),1)*[1,-2,1],-1:1,m(i),m(i))/h(i)^2;
D([1,end]) = -D([2,end-1]);

function runMinimalExample
% --- 2D --- 
omega = [0 2 0 2];
m     = [17 19];
xc    = getCellCenteredGrid(omega,m);
uc    = 1e-2*randn(size(xc));
[mbS, mbdS, mbd2S] = feval(mfilename,uc,omega,m,'alpha',1,'matrixFree',0);
[mfS, mfdS, mfd2S] = feval(mfilename,uc,omega,m,'alpha',1,'matrixFree',1);
fprintf('Curvature regularizer : %f (mb) | %f (mf) \n' , mbS,mfS);
RE = norm(mbdS-mfdS)/norm(mbdS);
fprintf('rel.error dS (mb vs. mf) : %e \n' , RE);
figure(1)
subplot(1,2,1);
spy(mbd2S)
title(sprintf('%s-spy(d2S)-(2D)',mfilename));
% --- 3D --- 
omega = [0 2 0 3 0 2];
m     = [6 5 7];
xc    = getCellCenteredGrid(omega,m);
uc    = 1e-2*randn(size(xc));
[mbS, mbdS, mbd2S] = feval(mfilename,uc,omega,m,'alpha',1,'matrixFree',0);
[mfS, mfdS, mfd2S] = feval(mfilename,uc,omega,m,'alpha',1,'matrixFree',1);
fprintf('Curvature regularizerl : %f (mb) | %f (mf) \n' , mbS,mfS);
RE = norm(mbdS-mfdS)/norm(mbdS);
fprintf('rel.error dS (mb vs. mf) : %e \n' , RE);
figure(1)
subplot(1,2,2);
spy(mbd2S)
title(sprintf('%s-spy(d2S)-(3D)',mfilename));
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
