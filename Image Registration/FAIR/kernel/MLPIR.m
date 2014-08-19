%==============================================================================
% (c) Jan Modersitzki 2011/01/04, see FAIR.2011 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%==============================================================================
%
% function [wOpt,his] = MLPIR(MLdata,varargin)
%
% Multi-Level Parametric Image Registration
%
% minimizes J^h(w) = D^h(T(Y(w)),R) + S^h(w) for h=coarse:fine
% uses PIR [=@GaussNewton]  for optimization
%               
% Input:
%   MLdata      struct of coarse to fine representations of the data, 
%               see getMultilevel.m
%  varagin     optinonal parameters, see default parameter
% Output:
%  wOpt         optimal parameter for the parametric part
%  MLhis        iteration history
%
% for level=minLevel:maxLevel,
%   get data(level)
%   if level>minLevel, w0 = wOpt; end;
%   get wOpt running PIR  using w0 as starting guess
% end%For
%
%==============================================================================

function [wOpt,his] = MLPIR(MLdata,varargin)

if nargin == 0, help(mfilename); E9_3Dbrain_MLPIR_sparse; return; end;

% setup default parameter for parametric pre-registration
PIRobj      = @PIRobjFctn;  % objective function for PIR
PIRopt      = @GaussNewton; % optimizer
PIRLS       = @Armijo;      % linesearch scheme
w0          = trafo('w0');  % starting guess
w0=[0 30 30]'
wStop       = w0;           % global stopping for PIR
wRef        = w0;           % regularization: (w-wRef)'*M*(w-wRef)
M           = [];           %
beta        = 0;            % regularization: H -> H + beta*I
maxIter     = 100;           % maximum number of iterations for PIR   
solver      = '';           % solver for PIR
getGrid     = @getCellCenteredGrid;
% setup additional default parameter
pause       = 0;            % flag for pauses
plots       = 1;            % flag for plots
plotIter    = 0;            % flag for output of iteration history each level
plotMLiter  = 0;            % flag for output of summarized iteration history
dimstr      = @(m) sprintf('m = [%s]',sprintf(' %d',m));
[MLdata,minLevel,maxLevel] = getMultilevel(MLdata);

for k=1:2:length(varargin),   % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
% build parameter cell for PIR
PIRpara = setOptPara(...
  'optimizer',  PIRopt,...
  'lineSearch', PIRLS,...
  'maxIter',    maxIter,...
  'Plots',      @FAIRplots,...
  'yStop',      wStop,...
  'solver',     solver,...
  'tolG',       1e-2,...
  varargin{:});

% FAIRmessage = @(str) fprintf('%% %s  [ %s ]  % s\n',...
%   char(ones(1,10)*'='),str,char(ones(1,60-length(str))*'='));
reportStatus
FAIRmessage(sprintf('%s, minLevel=%d:maxLevel=%d',...
  'MultiLevel Parametric Image Registration',minLevel,maxLevel));
omega = MLdata{end}.omega;
wOpt  = w0;
his   = [];
tic;

% -- for loop over all levels ---------------------------------------------
for level=minLevel:maxLevel;

  FAIRmessage(sprintf('%s: level %d from %d to %d, %s',...
    mfilename,level,minLevel,maxLevel,dimstr(MLdata{level}.m)));

  
  % get data for current level, compute interpolation coefficients
  m     = MLdata{level}.m; 
  [T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
  % update transformation
  trafo('set','omega',omega,'m',m);

  % initialize plots
  FAIRplots('reset','mode','PIR-multi level','fig',level,'plots',plots);
  FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

  % ----- call PIR ------------------------------------
  xc   = getGrid(omega,m); 
  Rc   = inter(R,omega,center(xc,m));
  fctn = @(wc) PIRobj(T,Rc,omega,m,beta,M,wRef,xc,wc);
  if level == minLevel, 
    fctn([]);   % report status
  else
    w0 = wOpt;  % update starting guess
  end; 
  [wOpt,hisPIR] = optim(PIRpara{:},'objFctn',fctn,'yc',w0,'yStop',wStop);
  % ----- end PIR --------------------------------------

  if plotIter,                                                
    plotIterationHistory(hisPIR,'J',[1,2,5],'fig',20+level);  
  end;                                                        
  
  
  % update iteration history
  if level == minLevel,
    his.str = hisPIR.str;
    his.his = hisPIR.his;
  else
    his.his = [his.his;hisPIR.his];
  end;
  doPause(pause)
  
end;%for level
% -- for loop over all levels ---------------------------------------------
his.time = toc;

% if plotMLiter,
%   plotMLIterationHistory(his,'fig',30);
% end;
if isempty(wStop), wStop = w0; end
his.reduction = fctn(wOpt)/fctn(wStop);
J = find(his.his(:,1)==-1); 
his.iter(minLevel:maxLevel) = his.his([J(2:end)-1;size(his.his,1)],1)';

FAIRmessage([mfilename,' : done !']);

%==============================================================================

function doPause(p)
if strcmp(p,'on'), 
  pause; 
elseif p>0, 
  pause(p); 
end;

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