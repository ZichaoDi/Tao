% (c) Jan Modersitzki 2011/04/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% function [fig,ph,th] = plotIterationHistory(his,varargin)
% 
% plots iteration history
% 
% Example: plotIterationHistory(his,'J',[1,2,5],'fig',2);
% Input:
%   his		    struct wigth headlines and list
%  varargin		optional parameters, see below
% Output:
%   fig         figure handle
%   ph          plot handle
%   th          text handle

function [fig,ph,th] = plotIterationHistory(his,varargin)

if nargin==0
    runMinimalExample; return;
end

fig     = [];
col     = 'krgbcmykrgbcmy';
figname = mfilename;
J       = [];
for k=1:2:length(varargin),
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(J), J = 1:length(his.str);  end;
str = his.str(J); his = his.his(:,J);
dum = abs(max(his,[],1));
dum = dum + (dum==0);
if ~isempty(fig),
  fig = figure(fig); clf; 
else  
  fig = figure;      clf; 
end;

set(fig,'numbertitle','off','name',sprintf('[FAIR:%d] %s',fig,figname));
for j=2:size(his,2),
  ph(j-1) = plot(his(2:end,1),his(2:end,j)/dum(j),'color',col(j-1)); hold on;
end;
set(ph,'linewidth',2);
th(1) = title(str{1});
th(2) = xlabel('iter');
th(3) = legend(ph,str{2:end});

function runMinimalExample
clear
close all
% setup data, interpolation, transformation, distance, 
% extract data from a particular level amd compute coefficients for interpolation
setupHandData;
inter('reset','inter','splineInter','regularizer','moments','theta',1e-1);
level = 4; omega = MLdata{level}.omega; m = MLdata{level}.m;
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega,'out',0);
distance('reset','distance','SSD');
trafo('reset','trafo','rotation2D','c',(omega(2:2:end)-omega(1:2:end))'/2); 
w0 = trafo('w0'); beta = 0; M =[]; wRef = []; % disable regularization
xc = getCellCenteredGrid(omega,m); 
Rc = inter(R,omega,xc);
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc);
% optimize
[~,his] = GaussNewton(fctn,w0,'plots',0);

plotIterationHistory(his);
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