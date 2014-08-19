% (c) Jan Modersitzki 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Test transformation schemes

toolbox = fileparts(which('trafo'));
testStart

debit = {
  'contents.m'
  'trafo.m'
  'affine2D.m'
  'affine2Dsparse.m'
  'affine3D.m'
  'affine3Dsparse.m'
  'matVecQw.m'
  'rigid2D.m'
  'rigid3D.m'
  'rotation2D.m'
  'splineTransformation2D.m'
  'splineTransformation2Dsparse.m'
  'splineTransformation3Dsparse.m'
  'translation2D.m'
  'translation3D.m'
  'tensorProdC.c'
  };


checkToolbox(toolbox,debit,'EDIT',0,'KEYBOARD',0,'RUN',1)

%% test syntax of trafo
trafo('clear');
trafo('disp');
trafo('reset','trafo','splineTransformation2D',...
  'p',[4,5],'omega',[1,1],'m',[6,40]);
trafo('disp');
[scheme,parameter] = trafo
pause(2);


%% test 2D cases

omega = [0,1,0,2]; m = [6,7]; xc = getCellCenteredGrid(omega,m);
center = (omega(2:2:end)+omega(1:2:end))'/2;
trafos = {
  'affine2D',{},...
  'affine2Dsparse',{},...
  'rigid2D',{},...
  'rotation2D',{'c',center},...
  'translation2D',{},...
  'splineTransformation2D',{'p',[4,5],'omega',omega,'m',m},...
  'splineTransformation2Dsparse',{'p',[4,5],'omega',omega,'m',m}
  };

for k=1:length(trafos)/2,
  fprintf('\n\n\ntest [ %s ] --- %d of %d\n',...
    trafos{2*k-1},k,length(trafos)/2)
  optn = trafos{2*k};
  trafo('reset','trafo',trafos{2*k-1},'debug','on',optn{:});
  trafo('disp');
  wc = trafo('w0');
  fctn = @(wc) trafo(wc,xc);
  checkDerivative(fctn,wc+randn(size(wc)));
  title(trafo)
  pause(2);
end;

%% test 3D cases
omega = [0,1,0,2,0,3]; m = [6,7,8]; xc = getCellCenteredGrid(omega,m);
trafos = {
  'affine3D',{},...
  'affine3Dsparse',{},...
  'rigid3D',{},...
  'splineTransformation3Dsparse',{'p',[4,5,6],'omega',omega,'m',m},...
  'translation3D',{}
  };

for k=1:length(trafos)/2,
  fprintf('\n\n\ntest [ %s ] --- %d of %d\n',...
    trafos{2*k-1},k,length(trafos)/2)
  optn = trafos{2*k};
  trafo('reset','trafo',trafos{2*k-1},'debug','on',optn{:});
  trafo('disp');
  wc = trafo('w0');
  fctn = @(wc) trafo(wc,xc);
  checkDerivative(fctn,wc+randn(size(wc)));
  title(trafo)
  pause(2);
end;

testEnd

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