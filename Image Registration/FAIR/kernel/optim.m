% ===============================================================================
% (c) Jan Modersitzki 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Main function for optimization, possible calls:
% 
% NPIRfctn = @(yc) NPIRobj(T,Rc,omega,m,yRef,yc); 
% NPIRpara = setOptPara(...
%  'optimizer', @GaussNewton,...
%  'lineSearch',@Armijo,...
%  'maxIter',   10,...
%  'Plots',     @FAIRplots,...
%  'yStop',     yStop,...
%  'solver',    'MATLAB'
% };
% yc = optim(NPIRpara{:},'objFctn',NPIRfctn,'yc',y0);
% run this file or see MLIR for an example.  
%==============================================================================

function [yc,His] = optim(varargin)

if nargin == 0, help(mfilename); minimalExample; return; end;

para = setOptPara(varargin{:});

[optimizer,para] = getValue(para,'optimizer');
[objFctn,para]   = getValue(para,'objFctn');
[yc,para]        = getValue(para,'yc');
[yc,His]         = optimizer(objFctn,yc,para{:});


function [value,optn] = getValue(optn,field)
j = find(~cellfun(@isempty,strfind(optn(1:2:end),field)),1,'last');
value = optn{2*j};
optn([2*j-1,2*j]) = [];

function minimalExample

% run minimnal example

% setup data and configure toolboxes building blocks
setupHandData;
viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
inter('reset','inter','splineInter','regularizer','moments','theta',1e-1);
level = 5; omega = MLdata{level}.omega; m = MLdata{level}.m;
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega,'out',0);
distance('reset','distance','SSD');
center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rotation2D','c',center); 
w0 = trafo('w0'); beta = 0; M =[]; wRef = []; % disable regularization

% initialize plots
FAIRplots('reset','mode','PIR-GN','fig',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

xc   = getCellCenteredGrid(omega,m); 
Rc   = inter(R,omega,xc);
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc);


% configure optimization parameter, 
PIRpara = setOptPara(...
 'optimizer', @GaussNewton,...
 'lineSearch',@Armijo,...
 'maxIter',   15,...
 'Plots',     @FAIRplots,...
 'yStop',     w0,...
 'solver',    'MATLAB');

fprintf('standard optimization parameters and updates\n')
disp(cell2struct(PIRpara(2:2:end),PIRpara(1:2:end),2))

wOpt = optim(PIRpara{:},'objFctn',fctn,'yc',w0,'yStop',w0);
fprintf('numerical optimizer: wOpt=%s\n',num2str(wOpt));

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