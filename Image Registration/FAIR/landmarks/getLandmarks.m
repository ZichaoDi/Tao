% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function [LM,fig] = getLandmarks(T,R,omega,m,varargin);
%
% Sets a number q of landmark in 2D template and reference images
%
% Input:
%   T,R         representation of template and reference
%   omega, m    representation of discretization for the domain
%   varargin    optional parameters
%
% Output:
%   LM 			q-by-4, the landmarks
%
% see also E5_Hands_TPS for an example

function [LM,fig] = getLandmarks(T,R,omega,m,varargin)

if nargin==0
  caller = dbstack;
  if any(strcmp({caller(:).name},'testLM')), return; end; % return in debug mode
  runMinimalExample; return;
end

LM = [];

% prepare visualization of T and R
omegaT = omega(1,:);
omegaR = omega(end,:);
xT  = getCellCenteredGrid(omegaT,m);
xR  = getCellCenteredGrid(omegaR,m);
Tc = inter(T,omegaT,xT,'inter','linearInter');
Rc = inter(R,omegaR,xR,'inter','linearInter');

fig = figure; clf; 
set(fig,'numbertitle','off','name',sprintf('[FAIR:%d]',fig));
subplot(1,2,1); viewImage(Tc,omegaT,m); hold on; title(inputname(1));
subplot(1,2,2); viewImage(Rc,omegaR,m); hold on; title(inputname(2));

% get a number of landmarks in T
t = getLM(Tc,omegaT,m,inputname(1),'r',fig,[1,2,1],inf);
if isempty(t), return;  end;

% get the same number of landmarks in R
r = getLM(Rc,omegaR,m,inputname(2),'g',fig,[1,2,2],size(t,1));
if isempty(r) || any(size(t) ~= size(r)), return;  end;
LM = [t,r];

function LM = getLM(Ic,omega,m,Iname,col,fig,sub,maxNumberLM);

fprintf('set landmarks in %s %s : % s\n',...
  Iname,'[left|middle|right]','[take|abort all|done]');
numberLM = 0;
while 1
  subplot(sub(1),sub(2),sub(3));
  [x,y,b] = ginput(1);
  switch b,
    case 1, % add new landmark
      numberLM = numberLM + 1;
      LM(numberLM,1:2) = [x,y];
      plot(x,y,[col,'.'],x,y,[col,'o'],'linewidth',2)
      text(x,y,int2str(numberLM),'color',col,'fontsize',20);
    case 2, % abort
      LM = [];
      return;
    case 3, % return
      return;
  end;
  % for the second image
  if numberLM == maxNumberLM, return; end;
end;

function runMinimalExample
setupHandData
getLandmarks(dataT,dataR,omega,m);
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