% ===============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%==============================================================================
%
% function para = setOptPara(varargin)
%
% sets default parameters for optimizers, see optim for an example
%==============================================================================

function para = setOptPara(varargin)

if nargin > 0 && isstruct(varargin{1}),
  para = varargin{1};
  varargin(1) = [];
else
  para = [];
  if nargout == 0,
    help(mfilename); 
    return;
  end;
end;
%varargin{:}

fields = {...
  'optimizer', '@GaussNewton',...
  'lineSearch','@Armijo',...
  'objFctn',[],...
  'solver',[],...
  'yc',[],...
  'yStop',[],...
  'JStop',[],...
  'vecNorm',@norm,... 
  'maxIter',10,...
  'tolJ',1e-4,...
  'tolY',1e-4,...
  'tolG',1e-2,...
  'LSMaxIter',10,...
  'LSreduction',1e-4,...
  'Plots',@(iter,para) [],...
  'plots',1
  }; 

% initialize para
for j=1:length(fields)/2
  if ~isfield(para,fields{2*j-1}),
    para = setfield(para,fields{2*j-1},fields{2*j});
  end;
end;
fields = {fields{1:2:end}};
para = orderfields(para,fields);

% update para
OK = 1;
for j=1:length(varargin)/2
  if ~isfield(para,varargin{2*j-1}),
    fprintf('warning: invalid parameter [%s]\n',varargin{2*j-1});
    OK = 0;
  else    
    para = setfield(para,varargin{2*j-1},varargin{2*j});
  end;
end;

para = reshape([fieldnames(para),struct2cell(para)]',1,[]);

return

% parameter initialization -----------------------------------------------
vecNorm      = @norm;           % norm to be used for dJ and dy    
solver       = [];              % linear solver 
yStop        = [];              % used for stopping in multi-level framework
Jstop        = [];              % 
Plots        = @(iter,para) []; % for plots;
for k=1:2:length(varargin),     % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
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