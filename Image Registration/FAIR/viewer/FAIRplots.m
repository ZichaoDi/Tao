% (c) Jan Modersitzki 2011/04/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% function FAIRplots(task,varargin);
% 
% FAIRplots: generates FAIR plots, see E6_FAIRplots for an example
%
% Create the following 2-by-3 plot:
% _________________________________________________________________
% | R0                 | T0                 |  Tk                 |
% | titleR0            | titleT0            |  title(Tk)          |
% | showImage(R0)      | showImage(T0)      |  showImage(Tk)      |
% |                    |                    |                     |
% |_______________________________________________________________|
% |                    | D0                 |  Dk                 |
% | titleGrid          | titleT0            |  title(Tk)          |
% | showImage(T0) G0   | showDifference(D0) |  showDifference(Tk) |
% | and showGrid  Gk   |                    |                     |
% |_______________________________________________________________|
%
% The function uses a persitent struct plotOptn controling its behaviour
% and perform the following tasks:
%
%------------------------------------------------------------------------------------
% - FAIRplots('clear')
%   clears all persistent variables and returns
%------------------------------------------------------------------------------------
% - FAIRplots('reset',varargin), FAIRplots('set',varargin)
%   initializes various variables and plot-handles (see details below)
%   overwrites defaults by the varargin list
%   note: no plots are shown using set
%   Example: FAIRplots('reset','mode','NPIR','fig',level);
%
%------------------------------------------------------------------------------------
% - FAIRplots('init',struct('Tc',dataT,'Rc',dataR,'omega',omega,'m',m));
%   initializes the figure 'fig' at position 'position, 
%     R(xc)            -> (2,3,1) title: Rname
%     T(xc)            -> (2,3,2) title: 'T(xc)'
%     T(xc)            -> (2,3,4) title: 'T(xc)'
%     T(xc)-R(xc)      -> (2,3,5) title: '|T(xc)-R|'
%------------------------------------------------------------------------------------
%   para  = struct('Tc',Tstop,'Rc',Rc,'omega',omega,'m',m,'yc',yStop,'Jc',Dc(Tstop));
% - FAIRplots('stop',para); T = Tstop 
%     T(yStop)         -> (2,3,3) title: 'Tstop'
%     T(yStop)-R(xc)   -> (2,3,6) title: '|Tstop-R|=100%'
%     yStop            -> (2,3,4) (update grid)
%------------------------------------------------------------------------------------
%  para.Tc = T0; para.yc = Y0; para.Jstop = para.Jc; para.Jc = Dc(T0);
% - FAIRplots('start',para); 
%     T(y0)            -> (2,3,2) title: Tname
%     T(y0)-R(xc)      -> (2,3,5) title: Dname
%     y0               -> (2,3,4) (update grid)
%------------------------------------------------------------------------------------
%   para.Tc = Tk; para.yc = yk; para.Jc = Dc(Tk); 
%   para.normdY = norm(yk-Y0)/norm(Ystop);
% - FAIRplots(k,para); (note: k is an integer)
%   T(yc)            -> (2,3,3) title: Tanme
%   T(yc)-R(xc)      -> (2,3,6) title: Dname
%   yc               -> (2,3,4) (replace)
%
%------------------------------------------------------------------------------------
%
%  Details on initialization (default values in []):
%
%   the figure is controlled by FIG[=[]] (handle), position[=screen dependent] (positioning),
%   subplots are controllable via handles [TRG]?HANDLE
%   PLOTS[=1] is used to disable the functionality (set plots=0)
%   MODE[=name of calling function] is used to control the default figurename and reference title:
%     mode==PIR*:   
%       figname: PIR*: inter/distance/trafo, dimension/size
%       Rname:   R, m=value of m, length(w)=number of parameters
%     mode==NPIR*:   
%       figname: NPIR*: inter/distance/regularizer, dimension/size
%       Rname:   R, m=value of m, alpha==value of alpha
%     otherwise
%       figname: mode, dimension/size
%       Rname:   R, R, m=value of m
%
%    omega and m are used at various places (computation and output)
%
% The final plots and their entiteling  are performed using reconfigurable function
%  
%    Function(arguments)  default
%    Tshow(T,omega,m)     @viewImage
%    Tname(k)             @(iter) sprintf('T(%d)',iter)
%                         
%    Rshow(R,omega,m)     @viewImage
%    Rname(m)             context dependent, 
%                         e.g. @(m)sprintf('R, %s, length(w)=%d',dimstr(m),wLen);
%                         
%    Dshow(T,R,omega,m)   @(T,R,omega,m) viewImage(255-abs(T-R),omega,m)
%    Dname(k,Jc,Jstop)    @(k,Jc,Jstop) ...
%                         sprintf('|J(%d)/Jstop|=%s%%',k,num2str(abs(100*Jc/Jstop)))
%    
%    Gshow(yc,omega,m)    context dependent, e.g. 
%                         @(yc,omega,m)  plotGrid(yc,omega,m,'spacing',ceil(m/32))
%    Gname(normdY)        @(normdY) sprintf('T(xc), |dY|= %s',num2str(normdY)) 
%    
% run('FAIRplots') or  - see also E6_FAIRplots for an example.

function FAIRplots(task,varargin)

if nargin==0
    runMinimalExample; return;
end

persistent plotOptn
task  = lower(task);

Error  = @(str) error(sprintf('[%s-task=%s] %s\n',mfilename,task,str));
dimstr = @(m) sprintf('[%s]',sprintf(' %d',m));

% ----- CLEAR plotOptn ----------------------------------------------
if strcmp(task,'clear') 
  if nargin == 1,
    clear plotOptions;
    disp('cleared plotOptn')
  else
    Error('nargin>0')
  end;
  return;
end;

% ----- INITIALIZE --------------------------------------------------
if strcmp(task,'reset') | strcmp(task,'set')
  % RESET plotOptn
  if strcmp(task,'reset'),
    plotOptn = [];
  end;

  % INITIALIZE plotOptn
  fields = {
    'fig'
    'plots'
    'figname'
    'position'
    'mode'
    'omega'
    'm'
    'T0handle'
    'Tkhandle'
    'Tshow'
    'Tname'
    'R0handle'
    'Rshow'
    'Rname'
    'D0handle'
    'Dkhandle'
    'Dshow'
    'Dname'
    'G0handle'
    'Gkhandle'
    'Gshow'
    'Gname' 
    'scale'
    };F(20) = struct('cdata',[],'colormap',[]);


  for j=1:length(fields)
    if ~isfield(plotOptn,fields{j}),
      plotOptn = setfield(plotOptn,fields{j},[]);
    end;
  end;

  % prepare the figname
  J = find(strcmp(varargin,'mode'));
  if isempty(J),
    mode = dbstack;               
    mode = mode(min(length(mode),2)).name;
  else
    mode = varargin{min(J)+1};
  end;

  % set defaults
  defaults = {
    'mode',     mode,...
    'plots',    1,...
    'position', [],...
    'Rshow',    @(R,omega,m) viewImage(R,omega,m,'scale',1),...     % 
    'Tshow',    @(T,omega,m) viewImage(T,omega,m,'scale',1),...%
    'Tname',    @(iter) sprintf('T(%d)',iter),...
    'Dshow',    @(T,R,omega,m) viewImage(255-abs(T-R),omega,m,'scale',1),...%
    'Dname',    @(j,Jc,J0) ...
        sprintf('|J(%d)/J_{0}|=%s%%',j,num2str(abs(100*Jc/J0))),...
    'Gshow',    @(yc,omega,m)  plotGrid(yc,omega,m,'spacing',ceil(m/32)),...
    'Gname',    @(normdY) sprintf('T(xc), |dY|= %s',num2str(normdY)) 
    };
 
  for j=1:2:length(defaults),
    if ~isfield(plotOptn,defaults{j}),
      plotOptn.(defaults{j}) = defaults{j+1};
    elseif isempty(plotOptn.(defaults{j})),
      plotOptn.(defaults{j}) = defaults{j+1};
    end;
  end;

  % collect information about current status of the toolbox
  try, wLen = length(trafo('w0'));      catch, wLen = nan;           end;
  try, interStr = inter;       catch, interStr = 'inter:void';       end;
  try, trafoStr = trafo;       catch, trafoStr = 'trafo:void';       end;
  try, distStr  = distance;    catch, distStr  = 'distance:void';    end;
  try, regStr   = regularizer; catch, regStr   = 'regularizer:void'; end;

  % if necessary TUNE figname
  if ~any(strcmp(varargin,'figname')),
  
    if length(plotOptn.mode) > 2 & strcmp(lower(plotOptn.mode(1:3)),'pir'),
      str = sprintf('%s: %s/%s/%s, ',plotOptn.mode,interStr,distStr,trafoStr);
    elseif length(plotOptn.mode) > 3 & strcmp(lower(plotOptn.mode(1:4)),'npir'),
      str = sprintf('%s: %s/%s/%s, ',plotOptn.mode,interStr,distStr,regStr);
    else
      str = sprintf('%s: ',mode);
    end;
    plotOptn.figname = ...
      @(omega,m) sprintf('%s%dD, m=%s',str,length(omega)/2,dimstr(m));;
  end;

  % if necessary TUNE Rname
  if ~any(strcmp(varargin,'Rname')),    
    if length(plotOptn.mode) > 2 & strcmp(lower(plotOptn.mode(1:3)),'pir'),
      plotOptn.Rname = @(m)sprintf('R, %s, length(w)=%d',dimstr(m),wLen);
    elseif length(plotOptn.mode) > 3 & strcmp(lower(plotOptn.mode(1:4)),'npir'),
      plotOptn.Rname = @(m) sprintf('R, %s, \\alpha=%s',...
        dimstr(m),num2str(regularizer('get','alpha')));
    else
      plotOptn.Rname = @(m) sprintf('R, %s',dimstr(m));
    end;
  end;
  
  % ----- OVERWRITE DEFAULTS ------------------------------------------
  for k=1:2:length(varargin), % overwrites defaults
    if ~isfield(plotOptn,varargin{k}),
      Error(sprintf('field %s not in use',varargin{k}));
    end;
    plotOptn.(varargin{k}) = varargin{k+1};
  end;
  % -------------------------------------------------------------------
  return;
end;
% -------------------------------------------------------------------


if ~plotOptn.plots, return; end;

fig = plotOptn.fig;

% extract the parameters for plots
if nargin>1, 
  T      = getField(varargin{1},'Tc'); 
  R      = getField(varargin{1},'Rc'); 
  Rsub=getField(varargin{1},'Rsub');
  omega  = getField(varargin{1},'omega');
  omegat  = getField(varargin{1},'omegat');
  m      = getField(varargin{1},'m'); 
  mt  = getField(varargin{1},'mt');
  yc     = getField(varargin{1},'yc');
  normdY = getField(varargin{1},'normdY');
  Jc     = getField(varargin{1},'Jc');
  if(isempty(mt)); mt=m;end
  if(isempty(omegat)); omegat=omega;end
end;

switch task,

  case 'clear', return; 
  case 'reset', return; 
  case 'set',   return; 
    
  case 'init',
    % activate figure, set colordef and figname
    
    if isempty(fig), fig = figure; else fig=figure(fig); end;
    plotOptn.fig = fig;
    % extract variables
    xc = getCellCenteredGrid(omega,m);
    xt = getCellCenteredGrid(omegat,mt);
%     T  = inter(reshape(T,m),omega,xc);
%     R  = inter(reshape(R,m),omega,xc);
    T  = inter(T,omegat,xt);
    R  = inter(R,omega,xc);

    if ~isempty(plotOptn.position),
      pos = plotOptn.position;
    else
      pos = FAIRposition('fig',fig);
    end;
    FAIRfigure(fig,...
      'figname',plotOptn.figname(omega,m),...
      'position',pos);

    if size(omega,2) == 6, % disable Gshow
      plotOptn.Gshow = @(yc,omegat,mt) [];
    end;

    % plot
    % _____________________________
    % | R       | T     |         |
    % |_________|_______|_________|
    % | T +grid | T -R  |         |
    % |_________|_______|_________|
   
    figure(fig); clf
      subplot(2,3,1);
        plotOptn.R0handle = plotOptn.Rshow(R,omega,m);      
        title(plotOptn.Rname(m));
      subplot(2,3,2);
        plotOptn.T0handle = plotOptn.Tshow(T,omegat,mt);  
        title('T(xc)');
      subplot(2,3,4); 
        plotOptn.G0handle = plotOptn.Tshow(T,omegat,mt);      
        plotOptn.Gkhandle = [];      
        title('T(xc)&grid'); hold on;
      subplot(2,3,5);
%%%============ plot based on optimizing on subimage 
      %R=reshape(R,m);
      %R=R(omegat(1)+1:omegat(2),omegat(3)+1:omegat(4));
      %plotOptn.D0handle = plotOptn.Dshow(T,R(:),omegat,mt); 
      Tnew=zeros(m);
      Tnew(omegat(1)+1:omegat(2),omegat(3)+1:omegat(4))=reshape(T,mt);
      plotOptn.D0handle = plotOptn.Dshow(Tnew(:),R(:),omega,m); 
%%%==========================================================
%          plotOptn.D0handle = plotOptn.Dshow(T,R,omegat,mt);    
        plotOptn.Dkhandle = [];      
        title('|T(xc)-R|');
    %pause(1/100); 
    drawnow;

  case 'stop',
    % plot
    % _____________________________
    % |         |       | Tstop   |
    % |_________|_______|_________|
    % |   +grid |       | Tstop-R |
    % |_________|_______|_________|
    figure(fig);
      subplot(2,3,3); 
        plotOptn.Tkhandle = plotOptn.Tshow(T,omegat,mt);     
        title('T^{stop}');
      subplot(2,3,6); 
            Tnew=zeros(m);
      Tnew(omegat(1)+1:omegat(2),omegat(3)+1:omegat(4))=reshape(T,mt);
      plotOptn.Dkhandle = plotOptn.Dshow(Tnew(:),R(:),omega,m); 
       % plotOptn.Dkhandle = plotOptn.Dshow(T,R,omegat,mt);
        title('|T^{stop}-R|, J^{0}=100%');
      subplot(2,3,4); 
        set(plotOptn.Gkhandle,'visible','off');
        plotOptn.Gkhandle = plotOptn.Gshow(yc,omegat,mt);
    %pause(1/100); 
    drawnow;
    plotOptn.Jstop = Jc;
    
  case 'start',
    % plot
    % _____________________________
    % |         | T0    |         |
    % |_________|_______|_________|
    % |   +grid | T0-R  |         |
    % |_________|_______|_________|
    figure(fig);
      subplot(2,3,2); 
        plotOptn.T0handle = plotOptn.Tshow(T,omegat,mt);    
        title(plotOptn.Tname(0));
      subplot(2,3,5); 
            Tnew=zeros(m);
      Tnew(omegat(1)+1:omegat(2),omegat(3)+1:omegat(4))=reshape(T,mt);
      plotOptn.D0handle = plotOptn.Dshow(Tnew(:),R(:),omega,m); 
      %  plotOptn.D0handle = plotOptn.Dshow(T,R,omegat,mt);   
        title(plotOptn.Dname(0,Jc,plotOptn.Jstop));
     subplot(2,3,4); 
       set(plotOptn.Gkhandle,'visible','off');
       plotOptn.Gkhandle = plotOptn.Gshow(yc,omegat,mt);
     drawnow;

  otherwise, 
    if ~isnumeric(task),
      warning(['don''t no how to deal task <',task,'>!']);
      return;
    end;
    
    % plot
    % _____________________________
    % |         |       | Tc      |
    % |_________|_______|_________|
    % |   +grid |       | Tc-R    |
    % |_________|_______|_________|    
    figure(fig)
      subplot(2,3,4); 
        set(plotOptn.Gkhandle,'visible','off');
        plotOptn.Gkhandle = plotOptn.Gshow(yc,omegat,mt); 
        title(plotOptn.Gname(normdY));
      subplot(2,3,3); 
        plotOptn.Tkhandle = plotOptn.Tshow(T,omegat,mt);    
        title(plotOptn.Tname(task));
      subplot(2,3,6); 
      Tnew=zeros(m);
      Tnew(omegat(1)+1:omegat(2),omegat(3)+1:omegat(4))=reshape(T,mt);
      plotOptn.Dkhandle = plotOptn.Dshow(Tnew(:),R(:),omega,m); 
%         plotOptn.Dkhandle = plotOptn.Dshow(T,R,omegat,mt);      
        title(plotOptn.Dname(task,Jc,plotOptn.Jstop));
    drawnow; 
 
end;
% figure(112);viewImage(T,omegat,mt);figure(113);viewImage(R,omegat,mt);

function v = getField(s,field);
if isfield(s,field), v = s.(field); else v = []; end;

function s = setDefault(s,varargin);
for j=1:2:length(varargin),
  if ~isfield(s,varargin{j}),
    s.(varargin{j}) = varargin{j+1};
  elseif isempty(s.(varargin{j})),
    s.(varargin{j}) = varargin{j+1};
  end;
end;

function runMinimalExample
setupHandData
inter('reset','inter','splineInter','regularizer','moments','theta',1e0);
level = 4; omega = MLdata{level}.omega; m = MLdata{level}.m;
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega(1,:),'out',0);
x0  = getCellCenteredGrid(omega(1,:),m);
x1  = rotation2D( 25*pi/180,x0,'c',(omega(2:2:end)-omega(1:2:end))'/2);
x2  = rotation2D(-25*pi/180,x0,'c',(omega(2:2:end)-omega(1:2:end))'/2);  
Rc  = inter(R,omega(end,:),x0);
T0  = inter(T,omega(end,:),x0);
T1  = inter(T,omega(end,:),x1);
T2  = inter(T,omega(end,:),x2);

FAIRfigure(2);
subplot(1,2,1); viewImage2Dsc(T0,omega,m,'title','template','colormap','gray')
subplot(1,2,2); viewImage2Dsc(Rc,omega,m,'title','reference','colormap','gray')

FAIRplots('reset','mode','testing-PIR','fig',10,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m));

para = struct('Tc',T0,'Rc',Rc,'omega',omega,'m',m,'yc',x0,'Jc',100);
FAIRplots('stop',para); 
%pause(3);

para = struct('Tc',T1,'Rc',Rc,'omega',omega,'m',m,'yc',x1,'Jc',14.3);
FAIRplots('start',para); 
%pause(3);

para = struct('Tc',T2,'Rc',Rc,'omega',omega,'m',m,'yc',x2,'Jc',0.3);
FAIRplots(14,para);
help(mfilename)

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