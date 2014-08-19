% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function [Y,LM] = LMreg(type,LM,X,varargin)
%
% Landmark based Image Registration
%
% Input:
%   type 		in {'linear','quadratic','TPS'} denotes the transformation model
%   LM          the location of landmarks m-by-2*d, d = dimension
%               t = LM(:,1:d), LM in template, r = LM(:,d+1:2*d) LM in reference
%   X           grid where the LM solution is evaluated
%   varargin    optional arguments, like theta for TPS, see below
%
% Output:
%    Y 			transformed grid, LM(:,5)
%    LM         updated landmarks
%
% see also E5_Hands_TPS for an example

function [Y,LM] = LMreg(type,LM,X,varargin)

if nargin==0
    runMinimalExample; return;
end

theta = 0;
for k=1:2:length(varargin), % overwrites default parameter 
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% reshape grid and extract landmarks for T and R
dim = size(LM,2)/2;
t   = LM(:,1:dim);
r   = LM(:,dim+1:2*dim);
X   = reshape(X,[],dim);
Y   = zeros(size(X));

switch type,
  case 'linear',
    % linear transformation model
    Q  = [r,ones(size(LM,1),1)];
    w  = (Q'*Q)\(Q'*t);
    for j = 1:dim,
      Y(:,j)        = [X,ones(size(X,1),1)]*w(:,j);
      LM(:,2*dim+j) = Q*w(:,j);      
    end;
  case 'quadratic',
    % quadratic 2D transformation model
    Q = [ones(size(LM,1),1),LM(:,[3:4]),LM(:,3).^2,LM(:,4).^2,LM(:,3).*LM(:,4)];
    w = (Q'*Q)\(Q'*LM(:,1:2));
    m = @(w,X) (w(1)+w(2)*X(:,1)+w(3)*X(:,2)...
      +w(4)*X(:,1).^2+w(5)*X(:,2).^2+w(6)*X(:,1).*X(:,2)); 
    for j = 1:dim,
      Y(:,j)        = [ones(size(X,1),1),X,X(:,1).^2,X(:,2).^2,X(:,1).*X(:,2)]*w(:,j);
      LM(:,2*dim+j) = Q*w(:,j);      
    end;
  case 'TPS',
    % Thin-Plate-Spline transformation model, 
    w = getTPScoefficients(LM(:,1:4),'theta',theta);
    [Y,yLM] = evalTPS(LM,w,X);
    LM(:,[5,6]) = reshape(yLM,[],2);
    
end;

if (nargout < 2) || (dim ~= 2), return; end;

% approximate the inverse transformation on the landmarks
Y = reshape(Y,[],2);
LM(:,7) = griddata(Y(:,1),Y(:,2),X(:,1),LM(:,1),LM(:,2));
LM(:,8) = griddata(Y(:,1),Y(:,2),X(:,2),LM(:,1),LM(:,2));
Y = Y(:);

function runMinimalExample
setupHandData
omegaT = omega(1,:);
omegaR = omega(end,:);
xT = getCellCenteredGrid(omegaT,m);
xR = getCellCenteredGrid(omegaR,m);
Tc = inter(dataT,omegaT,xT);
Rc = inter(dataR,omegaR,xR);

theta = 10;
[yc,LM] = LMreg('TPS',LM(:,1:4),xR,'theta',theta);
TLM = inter(dataT,omegaT,yc);

  
FAIRfigure(1,'figname',mfilename); clf;

subplot(1,3,1); 
  viewImage(Tc,omegaT,m); hold on; axis off;
  ph = plotLM(LM(:,1:2),'numbering','on','color','r');
  set(ph,'linewidth',2,'markersize',20);
  title(sprintf('%s','T&LM'),'fontsize',20);

subplot(1,3,2); 
  viewImage(Rc,omegaR,m); hold on;  axis off;
  ph = plotLM(LM(:,3:4),'numbering','on','color','g','marker','+');
  set(ph,'linewidth',2,'markersize',20);
  title(sprintf('%s','R&LM'),'fontsize',20);

subplot(1,3,3); 
  viewImage(TLM,omegaR,m); hold on;  axis off;
  ph = plotLM(LM(:,3:4),'numbering','off','color','g','marker','+');
  qh = plotLM(LM(:,7:8),'numbering','off','color','m','marker','x');
  rh = plot(LM(:,[3,7])',LM(:,[4,8])','m-','linewidth',3);
  set([ph;qh;rh],'linewidth',2,'markersize',20);
  title(sprintf('T(Y^{TPS},\\theta=%s)&LM',num2str(theta)),'fontsize',30);

help(mfilename);

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
