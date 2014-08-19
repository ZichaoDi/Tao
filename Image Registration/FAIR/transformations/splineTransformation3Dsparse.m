% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function [y,dy] = splineTransformation3Dsparse(w,x,varargin)
%
% computes y = dy(w) and returns dy = function_handle to 
% kron(I2,Q3,Q2,Q1), Q{i}(:,1) = spline(x(:,1));
% if no arguments are given, the parameters for the identity map are returned.
%
% required inputs are: p (number of spline coefficients), omega and m
% 
% see also transformations/contents.m, trafo.m and BigTutorialTrafo

function [y,dy] = splineTransformation3Dsparse(w,x,varargin)

% the persistent variable Q stores the matrix 
% Q(x) = kron( I_2 , Q3, Q2, Q1 );
% p is number of spline coefficients, m is size of grid, omega is domain
persistent Q p m omega

% p = []; m = []; omega = []; 
flag = 'y';

for k=1:2:length(varargin), % overwrites defaults
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if nargin == 0, 
  runMinimalExample;
  return;
else
  y = mfilename('fullfile'); 
  dy = zeros(3*prod(p),1);             % parameterization of identity
  if ischar(w),  Q  = []; w = []; end; % reset Q
  if isempty(w), return;          end; 
end;

OK = 0;
if isstruct(Q),
  if isfield(Q,'m'), mQ = Q.m; else mQ = 0*m; end;
  if isfield(Q,'p'), pQ = Q.p; else pQ = 0*p; end;
  OK = all(p == pQ & m == mQ);
end;

if ~OK,
  % it is assumed that x is a cell centered grid, extract xi1 and xi2
  x          = reshape(x,[m,3]);
  Q.omega    = omega;
  Q.m        = m;
  Q.dim      = 3;
  Q.p        = p;
  Q.Q1       = getQ1d(omega(1:2),m(1),p(1),x(:,1,1,1));
  Q.Q2       = getQ1d(omega(3:4),m(2),p(2),x(1,:,1,2));
  Q.Q3       = getQ1d(omega(5:6),m(3),p(3),x(1,1,:,3));
  Q.Q        = @(w) splineTransformation3Dsparse(w,x,'flag','Qw');
  Q.Qadjoint = @(w) splineTransformation3Dsparse(w,x,'flag','QTw');
  if nargout == 0, return; end;
end;
dy = Q;

if strcmp(flag,'Qw') || strcmp(flag,'y'),
  w  = reshape(w,[Q.p,Q.dim]);
  y  = tensorProd(w,Q.Q1,Q.Q2,Q.Q3);
  y  = y(:);
  if strcmp(flag,'y'), y  = x(:) + y(:); end;
else % flag is 'QTw'
  w  = reshape(w,[Q.m,Q.dim]);
  y  = tensorProd(w,Q.Q1',Q.Q2',Q.Q3');
  y  = y(:);
end

function Q = getQ1d(omega,m,p,xi)
Q  = zeros(m,p); xi = reshape(xi,[],1);
for j=1:p,
  cj=zeros(p,1); cj(j) = 1;
  Q(:,j) = splineInter(cj,omega,xi);
end;
Q = sparse(Q);


% the following functions all compute kron(Q3,Q2,Q1) * w which is
% computationally challenging
%
% while all functions give equivalent results, tensorProd is the fastest
% implementation
function Qw = tensorProd(w,Q1,Q2,Q3)
p = [size(Q1);size(Q2);size(Q3)];
Qw = zeros([p(:,1)',3]);

for k=1:3,
  w0 = w(:,:,:,k);
  w0 = reshape(Q1*reshape(w0,p(1,2),[]),p(1,1),p(2,2),p(3,2));
  w0 = permute(w0,[2,1,3]);
  w0 = reshape(Q2*reshape(w0,p(2,2),[]),p(2,1),p(1,1),p(3,2));
  w0 = ipermute(w0,[2,1,3]);
  w0 = permute(w0,[3,1,2]);
  w0 = reshape(Q3*reshape(w0,p(3,2),[]),p(3,1),p(1,1),p(2,1));
  w0 = ipermute(w0,[3,1,2]);
  Qw(:,:,:,k) = w0;
end;

function Qw = tensorProd2(w,Q1,Q2,Q3)
p = [size(Q1);size(Q2);size(Q3)];
Qw = zeros([p(:,1)',3]);

for k=1:3,
  for i1=1:size(Qw,1),
    for i2=1:size(Qw,2),
      for i3=1:size(Qw,3),
        for j1=1:size(w,1),
          for j2=1:size(w,2),
            for j3=1:size(w,3),
              Qw(i1,i2,i3,k) = Qw(i1,i2,i3,k) + Q1(i1,j1)*Q2(i2,j2)*Q3(i3,j3)*w(j1,j2,j3,k);
            end;
          end;
        end;
      end;
    end;
  end
end;

function Qw = tensorProd3(w,Q1,Q2,Q3) % C implemantation of tensorProd2
Qw = tensorProdC(w(:),Q1(:),Q2(:),Q3(:),[size(Q1,1) size(Q2,1) size(Q3,1)],size(w));

function runMinimalExample
help(mfilename);
fprintf('%s: minimal example\n',mfilename)

omega = [0,10,0,8,0,9]; m = [8,7,6]; p = [5,4,3];
w = zeros([p,3]);  w(3,3,3,1) = 0.05; w(3,4,3,2) = -0.01;
x = getCellCenteredGrid(omega,m);
y = feval(mfilename,w(:),x,'omega',omega,'m',m,'p',p,'Q',[]);
figure(2); clf;
plotGrid(x,omega,m,'color','r'); axis image; hold on;
plotGrid(y,omega,m,'color','b'); axis image; hold off;  

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