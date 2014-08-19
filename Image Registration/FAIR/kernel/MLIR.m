%==============================================================================
% (c) Jan Modersitzki 2011/01/04, see FAIR.2011 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%==============================================================================
%
% function [yc,wc,MLhis] = MLIR(MLdata,varargin)
%
% Multi-Level Image Registration
%
% minimizes J^h(y) = D^h(T(y),R) + S^h(y) for h=coarse:fine
% uses PIR [=@GaussNewton] and NPIR [=@GaussNewton]
%               
% Input:
%   MLdata    coarse to fine representations of the data, see getMultilevel.m
%   varagin     optinonal parameters, see default parameter
% Output:
%   yc          numerical optimizer
%   wc          optimal parameters for PIR
%   MLhis       iteration history
%
% for level=minLevel:maxLevel,
%   get data(level)
%   if level==minLevel, get wOpt running PIR; end;
%   use pre-registered Yref=trafo(wOpt,X) for regularization
%   if level==minLevel
%     y0 = Ystop; % use this as starting guess
%   else
%     get y0 by prologating Yopt from coarse to finer level
%   end;
%   get Yopt running NPIR using y0 as starting guess
% end%For
%
% see also E9_Hands_MLIR_SSD_mbElas for an example
%==============================================================================

function [yc,wc,MLhis] = MLIR(MLdata,varargin)

if nargin == 0, help(mfilename); E9_HNSP_MLIR_SSD_mfCurv; return; end;

% setup default parameter for parametric pre-registration
parametric  = 1;            % flag for parametric pre-registration
PIRopt      = @GaussNewton; % optimizer to be used for PIR
PIRLS       = @Armijo;      % optimizer to be used for PIR
PIRobj      = @PIRobjFctn;  % objective function for PIR
beta        = 0;            % regularization for Hessian in PIR
w0          = [];           % starting guess for PIR
wStop       = [];           % global stopping for PIR
wRef        = [];           % regularization: (w-wRef)'*M*(w-wRef)
M           = [];           %
maxIterPIR  = 50;           % maximum number of iterations for PIR
solverPIR   = '';           % solver for PIR

% setup default parameter for non-parametric registration
NPIRopt     = @GaussNewton; % optimizer to be used for NPIR
NPIRLS      = @Armijo;      % optimizer to be used for NPIR
NPIRobj     = @NPIRobjFctn; % objective function for NPIR
yStop       = [];           % global stopping for NPIR
yRef        = [];           % regularization: S(yc-yRef)
maxIterNPIR = 10;           % maximum number of iterations for NPIR
solverNPIR  = '';           % solver for NPIR

% setup additional default parameter
pause       = 0;        % flag for pauses
plots       = 1;        % flag for plots
plotIter    = 0;        % flag for output of iteration history each level
plotMLiter  = 1;        % flag for output of summarized iteration history

% initialize output
yc = []; wc = []; MLhis = [];

% get minLevel and maxLevel from MLdata
[MLdata,minLevel,maxLevel] = getMultilevel(MLdata);

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% build parameter cell for PIR
PIRpara = setOptPara(...
 'optimizer', PIRopt,...
 'lineSearch',PIRLS,...
 'maxIter',   maxIterPIR,...
 'Plots',     @FAIRplots,...
 'yStop',     wStop,...
 'solver',    solverPIR,...
 varargin{:});

% build parameter cell for PIR
NPIRpara = setOptPara(...
 'optimizer', NPIRopt,...
 'lineSearch',NPIRLS,...
 'maxIter',   maxIterNPIR,...
 'Plots',     @FAIRplots,...
 'yStop',     yStop,...
 'solver',    solverNPIR,...
 varargin{:});

% initialization
dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
omega   = MLdata{end}.omega;  % spatial domain
xc      = [];                 % current grid
grid    = regularizer('get','grid');
switch grid,
  case 'cell-centered', getGrid = @(m) getCellCenteredGrid(omega,m);
  case 'staggered',     getGrid = @(m) getStaggeredGrid(omega,m);
  case 'nodal',         getGrid = @(m) getNodalGrid(omega,m);
  otherwise, error('nyi');
end;

fprintf('\n\n');
fprintf('%s: MultiLevel Image Registration\n',mfilename)
fprintf('-- distance=%s, regularizer=%s, alpha=%s, trafo=%s\n',...
  distance,regularizer,num2str(regularizer('get','alpha')),trafo);

tic;
%--------------------------------------------------------------------------
for level=minLevel:maxLevel,

  FAIRmessage(sprintf('%s: level %d from %d to %d, %s',...
    mfilename,level,minLevel,maxLevel,dimstr(MLdata{level}.m)));

  % save(old grid), update(m,grid,data coefficients)
  xOld  = xc; 
  m     = MLdata{level}.m; 
  xc    = getGrid(m);
  [T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega);
  Rc    = inter(R,omega,center(xc,m));
    
  if level == minLevel && parametric, % parametric pre-registration
    % initialize plots
    FAIRplots('reset','mode','PIR','fig',minLevel-1,'plots',plots);
    FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

    % ----- call Parametric Image Registration ----------------------------
    xC = getGrid(m); 
    PIRfctn = @(wc) PIRobj(T,Rc,omega,m,beta,M,wRef,xC,wc); 
    PIRfctn([]); % report status of objective function
    if isempty(w0),    w0    = trafo('w0'); end;
    if isempty(wStop), wStop = w0;          end;
    
    [wc,his] = optim(PIRpara{:},'objFctn',PIRfctn,'yc',w0,'yStop',wStop);
    
    fprintf('wc = \n'); disp(wc');  doPause(pause);    
    % ----- end PIR --------------------------------------
    if plotIter,                              
      his.str{1} = 'iteration history for PIR'; 
      plotIterationHistory(his,'J',1:4,'fh',100+minLevel-1);
    end;                                      
  elseif level == minLevel, % no pre-registration
    wc = [];                                
  end;                                        
  
  % compute yRef = xc or yc = trafo(wc,xc) on the appropriate grid 
  if isempty(wc), % use yRef(x) = x
    yRef = xc;  
  else            % compute yn=y(wc,xNodal) and interpolate on current grid    
    yn   = trafo(wc,getNodalGrid(omega,m));   
    yRef = grid2grid(yn,m,'nodal',grid);  
  end;
    
  % compute starting guess y0
  if level == minLevel,
    y0 = yRef;  % the best known so far
  else
    % prolongate yc (coarse) y0 (current) 
    y0 = xc + mfPu(yc - xOld,omega,m/2);
  end;
  

  % ----- call NPIR -------------------------------------
  % S(y) = 0.5*alpha*hd*norm(B*(Y-Yref))^2
  % initialize plots for NPIR
  FAIRplots('reset','mode','NPIR','fig',level,'plots',plots);
  FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 
  NPIRfctn = @(yc) NPIRobj(T,Rc,omega,m,yRef,yc); 
  if level == minLevel, NPIRfctn([]);  end; % report status of objective function
  yStop = getGrid(m);
  [yc,his] = optim(NPIRpara{:},'objFctn',NPIRfctn,'yc',y0,'yStop',yStop);
  % ----- end NPIR --------------------------------------

  if plotIter,
    his.str{1} = sprintf('iteration history for NPIR, level=%d',level);
    plotIterationHistory(his,'J',1:4,'fh',100+level);
  end;

  % update iteration history
  if level == minLevel,
    MLhis.str = his.str;
    MLhis.his = his.his;
  else
    MLhis.his = [MLhis.his;his.his];
  end;
  doPause(pause);

%--------------------------------------------------------------------------
end;%For level
%--------------------------------------------------------------------------
MLhis.time = toc;
if plotMLiter,
  plotMLIterationHistory(MLhis,'fh',[]);
end;
if isempty(yStop), yStop = xc; end
MLhis.reduction = NPIRfctn(yc)/NPIRfctn(yStop);
J = find(MLhis.his(:,1)==-1); 
MLhis.iter(minLevel:maxLevel) = MLhis.his([J(2:end)-1;size(MLhis.his,1)],1)';

FAIRmessage([mfilename,' : done !']);

function doPause(p)
if strcmp(p,'on'), 
  FAIRpause; 
elseif p>0, 
  FAIRpause(p); 
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