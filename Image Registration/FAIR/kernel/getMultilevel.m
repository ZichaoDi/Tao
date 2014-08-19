%==============================================================================
% function [MLdata,minLevel,maxLevel,fig] = getMultilevel(IS,omega,m,varargin)
% (c) Jan Modersitzki 2009/07/30, see FAIR.2 and FAIRcopyright.m.
%
% compute multi-level representation of images IS = {T,R,MT,MR,...}
%
% Input:
%  IS        an image or a list of images
%  omega     describing the domain
%  m         number of discretization points
%  varargin  optional parameter, see below
%
% Output:
%   MLdata       struct containing multi-level representation
%    MLdata{j} =  {Tl,Rl,omega,ml}
%   minLevel     coarsest level
%   maxLevel     finest   level
%   fig          handle to the graphical output
%
% A continuous representation of an image serves as a starting point. Using a cell
% centered grid, the data is replace by interpolated values. The representation on
% a coarser level is obtained by averaging over adjacent cells (see filter for
% options). MLdata{level} is a structure containing the image(s), omega, and m where
% level runs 0:maxLevel, MLdata{level} = [] for level < minLevel.
%
%==============================================================================
% changes 2008-12-04
%  modified code to fit new omega definition
%==============================================================================

function [MLdata,minLevel,maxLevel,fig] = getMultilevel(IS,omega,m,varargin)

if nargin == 0,    % help and minimal example
    help(mfilename);
    Tdata = double(flipud(imread('hands-T.jpg'))');
    Rdata = double(flipud(imread('hands-R.jpg'))');
    omega = [0,20,0,25]; % specify physical domain
    m     = size(Tdata);
    % set view options and interpolation options and initialize viewer and interpolator
    viewOptn = {'viewImage','viewImage2D','colormap','gray(256)'};
    viewImage('reset',viewOptn{:});
    % create multilevel representation of the data
    MLdata = getMultilevel({Tdata,Rdata},omega,m,'fig',2);
    return;
end;

if nargin == 1,
    % simply return MLdata, minLevel, maxLevel
    MLdata = IS; maxLevel = length(MLdata); fig = 0;
    for minLevel=maxLevel:-1:1,
        if minLevel == 1 || isempty(MLdata{minLevel-1}), break; end;
    end;
    return;
end;

% start the work
fig      = 2;                  % figure number for output
dopause  = 0;                  % make pause for demonstrations
filter   = 'block';%'gaussian';%            % discrete smoothing kernel
names    = {'T','R','Q','Q1','Q2','Q3','Q4','Q5'};
restrictdim = ones(size(m));   % by default: restrict all dimensions
inter       =  @linearInter;

for k=1:2:length(varargin),    % overwrite defaults
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if ~iscell(IS), IS = {IS}; end;
lenIS = length(IS);                 % number of images to be handled
dim   = length(omega)/2;            % spacial dimension
maxLevel = ceil(log2(min(m)));      % finest level
minLevel = maxLevel-3;                  % minimal level size
msg = sprintf('%s(%d image(s), filter=%s, %dD data,level=%d:%d)',...
    mfilename,lenIS,filter,dim,minLevel,maxLevel);

fprintf('%s, figure=%d\n',msg,fig);

% do some plots
if fig,
    FAIRfigure(fig,'figname',msg);
end;

% some output to the console
fprintf('%s: [',mfilename);

omegak = @(k)omega(min(1+(k>1),size(omega,1)),:);

% start the loop from maxlevel to minlevel
for level = maxLevel:-1:minLevel
    if level == maxLevel,
        % replace data by the interpolant on x
        for k=1:lenIS,
            sizeData = size(IS{k});
            % note: data is either m1-by-....-by-md
            %                   or m1-by-....-by-md-by-nrVolume
            % but the 1d case is tricky
            
            if (dim>1) && (length(sizeData) > dim),
                nrVolume = sizeData(end),
                sizeData(end) = [];
                keyboard
            else
                nrVolume = 1;
            end;
            
            data  = reshape(IS{k},prod(sizeData),nrVolume);
            block = zeros(prod(m),nrVolume);
            xc = getCellCenteredGrid(omegak(k),m);
            for vol=1:size(data,2),
                block(:,vol) = inter(reshape(data(:,vol),sizeData),omegak(k),xc);
            end
            
            IS{k} = reshape(block,[m,nrVolume]);
        end;
        fprintf('%d',level);
    else
        % restrict the image to the coarser grid
        L = MLdata{level+1};
        for k=1:lenIS, % run over all images
            for j=1:dim, % run over all dimensions
                if restrictdim(j),
                    IS{k} = restrict(IS{k},dim,j,filter);
                end
            end;
        end;
        fprintf(',%d',level);
    end;
    
    % store the current level data in a struct
    L.m = size(IS{1}); L.m = L.m(1:dim); L.omega = omega;
    for k= 1:lenIS,
        L = setfield(L,names{k},IS{k});
    end;
    MLdata{level} = L;
    
    if fig, % do some plots
        str = @(k) sprintf('%s(level=%d), %s',...
            names{k},level,sprintf('m=[%s]',sprintf(' %d',L.m)));
        p0 = level-minLevel+1; dp = (maxLevel-minLevel+1);
        figure(fig);
        for k=1:lenIS,
            xc = getCellCenteredGrid(omegak(k),L.m);
            subplot(lenIS,dp,p0+(k-1)*dp);
            viewImage(inter(IS{k},omegak(k),xc),omegak(k),L.m);
            title(str(k));
        end;
        if dopause, pause; else drawnow; end;
    end;
end;
fprintf('] done\n');

%==============================================================================

function T = restrict(T,dim,j,filter)
J = [j,setdiff(1:dim,j)]; J = [J dim+1]; % bring dimension j to front
T = permute(T,J);
if rem(size(T,1),2), T(end+1,:,:,:) = T(end,:,:,:); end;
switch filter,
    case 'gaussian',
        Ta = zeros(size(T)+[2,zeros(1,ndims(T)-1)]);
        Ta(2:end-1,:,:,:) = T; Ta([1,end],:,:,:) = Ta([2,end-1],:,:,:);
        T =   0.125 * Ta(1:end-3,:,:,:) + 0.375 * Ta(2:end-2,:,:,:) ...
            + 0.375 * Ta(3:end-1,:,:,:) + 0.125 * Ta(4:end,:,:,:);
        T = T(1:2:end,:,:,:);
    case 'block'
        T = (T(1:2:end,:,:,:)+T(2:2:end,:,:,:))/2;
end;
T = ipermute(T,J); % bring j back home

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