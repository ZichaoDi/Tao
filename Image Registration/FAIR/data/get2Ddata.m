% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% Initializes 2D data to be used in FAIR.
%
% This file adds the following variables to outfile.mat:
%   dataT      data of template  image
%   dataR      data of reference image
%   omega      = [left,right,bottom,top]  domain description
%   m          size for finest data representation
%   MLdata     multilevel representation of data, see getMultilevel for details
%   viewOptn   parameterization of image viewer
%   intOptn  parameterization of image interpolater
%
%   Used in setup*data, see also data/contents.m

function get2Ddata(outfile,fileT,fileR,varargin)

if nargin == 0,
  % auto-test
  help(mfilename);
  return;
end;

if exist(outfile,'file'),
  fprintf('[%s] %s exists\n',mfilename,outfile);
  return
end;

omega  = [];
m      = [];

for k=1:2:length(varargin), % overwrite defaults  
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

txt = {
  'creates '
  ' - dataT:    data for templateimage, d-array of size p'
  ' - dataR:    data for reference image, d-array of size q'
  ' - omega:    coding domain = [omega(1),omega(2)]x...x[omega(2*d-1),omega(2*d)]'
  '             note: omega(1,:) for T' 
  '                   omega(end,:) for all others'
  ' - m:        size of the finest interpolation grid'
  ' - MLdata:   MLdata representation of the data'
  };
fprintf('%s\n',txt{:});

% do whatever needed to be done to get your data here
image = @(str) double(flipud(imread(str))'); % reads and converts

% load the original data, set domain, initial discretization, and grid
dataT = image(fileT); 
dataR = image(fileR);

if isempty(omega), omega = [0,size(dataT,1),0,size(dataT,2)];  end;
if isempty(m),     m     = size(dataR);                        end;
caller = dbstack;  caller = caller(min(2,length(caller))).name;

inter('reset','inter','linearInter');
intOptn = {'inter','linearInter'};
[viewer,viewOptn] = viewImage;

omegaI = @(i) omega(min(i,size(omega,1)),:);
xc = @(k) getCellCenteredGrid(omegaI(k),m);

viewData  = @(I,k) viewImage(inter(I,omegaI(k),xc(k)),omegaI(k),m);

FAIRfigure(1,'figname',caller); clf;
subplot(1,2,1); viewData(dataT,1); title('template');
subplot(1,2,2); viewData(dataR,2); title('reference');

MLdata = getMultilevel({dataT,dataR},omega,m,'fig',2);

% save to outfile
save(outfile,'dataT','dataR','omega','m','MLdata','viewOptn','intOptn');

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