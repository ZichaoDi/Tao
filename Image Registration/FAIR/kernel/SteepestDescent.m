% ===============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% ===============================================================================
%
% function [Yc,His] = SteepestDescent(fctn,Yc,varargin);
%
% Steepest Descent scheme with line search for minimizing J = fctn(Yc)
% 
% Input:
%   fctn         function handle
%   Yc           starting guess 
%   varargin     optional parameter, see below
%
% Output:
%   Yopt         numerical optimizer
%   His          iteration history
%
%==============================================================================

function [Yopt,His] = SteepestDescent(fctn,Yc,varargin)

if nargin ==0, help(mfilename); E6_Hands_PIR_SD; return; end;

% -- start parameter setup ----------------------------------------------
Ystop       = [];       % default reference parameter for stopping
maxIter     = 10;       % max number of iterations
tolJ        = 1e-3;     % stopping: relative variation of objective function
tolG        = 5e-2;     % stopping: relative variation of the gradient
tolY        = 5e-3;     % stopping: relative variation of the parameters
vecNorm     = @norm;
LSMaxIter   = 10;       % LineSearch: max number of line search steps
LSreduction = 1e-6;     % LineSearch: guaranteed reduction by line search
stepLength  = 1;        % scaling of gradient direction
Plots      = @(iter,para) []; % for plots;

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(Ystop), Ystop = Yc;          end;    % check Ystop
his    = zeros(maxIter+2,6);                    % summarizes iteration history

FAIRmessage = @(str) fprintf('%% %s  [ %s ]  % s\n',...
  char(ones(1,10)*'-'),str,char(ones(1,60-length(str))*'-'));
FAIRmessage(['JM 2008/08/06 : ',mfilename]);
% report objective function status and prepare plots
fctn([]);
fprintf('  %20s : %s\n','maxIter',num2str(maxIter));
fprintf('\n\n')
% -- end parameter setup   ----------------------------------------------


% -- start initial phase -------------------------------------------------
STOP = zeros(5,1);
[Jstop,para]  = fctn(Ystop);            % compute reference for stopping
Plots('stop',para);

iter = 0; Jold = 0; Yold = 0;           % set initial values
[Jc,para,dJ]  = fctn(Yc);               % compute current values
Plots('start',para);

his(1,:)  = [-1,Jstop,0,0,0,0]; 
his(2,:)  = [0,Jc,abs(Jc-Jstop),vecNorm(dJ),vecNorm(Yc-Yold),0]; 

% headlines for the output table
fprintf('%-4s %-12s %-12s %-12s %-12s %-4s\n',...
  'iter','J','Jold-Jc','|dJ|','|dy|','LS');
fprintf('%-4d %-12.2e %-12.2e %-12.2e %-12.2e %-4d\n',his(1,:));
% -- end initial phase   -------------------------------------------------

t = 1; % initial step
% -- start iteration phase -----------------------------------------------
while 1,
  tOld = t;
  iter    = iter + 1;
  
  % check stopping rules
  STOP(1) = abs(Jold-Jc)  <= tolJ*(1+abs(Jstop));
  STOP(2) = norm(Yc-Yold) <= tolY*(1+norm(Yc));
  STOP(3) = norm(dJ)      <= tolG*(1+abs(Jstop));
  STOP(4) = norm(dJ)      <= 1e3*eps;
  STOP(5) = (iter > maxIter);
  if all(STOP(1:3)) | any(STOP(4:5)), break;  end;
  
  % scaled steepest descent direction
  dY  = -stepLength*dJ';
  
  % perform line-search
  [t,Yt,LSiter] = Armijo(fctn,Yc,dY,Jc,dJ,....
    'LSMaxIter',LSMaxIter,'LSreduction',LSreduction);
  if (t == 0), break; end; % break if line-search fails
  
  if (t == 1),
    % everything is fine, be more optimistic
    stepLength = 2*stepLength;  
  else
    % adjust stepLength for the next step
    stepLength = stepLength*t;  
  end;
  
  % update   
  Yold = Yc; Jold = Jc;  Yc = Yt; Yold = Yc;
  [Jc,para,dJ] = fctn(Yc);
  para.normdY = vecNorm(Yc-Yold);
  Plots(iter,para);

  % store intermediate results, some output
  his(iter+2,:)  = [iter,Jc,abs(Jc-Jold),norm(dJ),para.normdY,LSiter]; 
  fprintf('%-4d %-12.2e %-12.2e %-12.2e %-12.2e %-4d\n',his(iter+1,:));
  
end;%while
% -- end iteration phase -------------------------------------------------

Yopt = Yc;
His.his = his(1:iter+1,:);
His.str = {'iter','J','Jold-Jc','|\nabla J|','|dY|','LS'};

fprintf('\n');
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e ]\n',STOP(1),...
  '(Jold-Jc)',(Jold-Jc),'tolJ*(1+|Jstop|)',tolJ*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e ]\n',STOP(2),...
  '|Yc-Yold|',norm(Yc-Yold),'tolY*(1+norm(Yc)) ',tolY*(1+norm(Yc)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e ]\n',STOP(3),...
  '|dJ|',norm(dJ),'tolG*(1+abs(Jstop)',tolG*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e ]\n',STOP(4),...
  'norm(dJ)',norm(dJ),'eps',1e3*eps);
fprintf('%d[ %-10s=%-16d  > %-25s=%-16d ]\n',STOP(5),...
  'iter',iter,'maxIter',maxIter);

FAIRmessage([mfilename,' : done !']);

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