% (c) Jan Modersitzki 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% FAIR-testing: the interpolation toolbox

% ----------------------------------------------------------------------------------------

toolbox = fileparts(which('inter'));
testStart;

debit = {
  'contents.m'
  'inter.m'
  'nnInter.m'
  'linearInter.m'
  'linearInterMex.m'
  'linearInterMexC.cpp'
  'linearInterSmooth.m'
  'linearInterSmoothMex.m'
  'linearInterSmoothMexC.cpp'
  'linearInterMatlab.m'
  'getSplineCoefficients.m'
  'motherSpline.m'
  'splineInter.m'
  'splineInterMex.m'
  'splineInterMexC.cpp'
  'cubicInter.m'                            
  'cubicInterMex.m'                         
  'cubicInterMexC.cpp'                      

  };

checkToolbox(toolbox,debit,'KEYBOARD',0,'EDIT',0,'MEX',1,'RUN',1)

%% test syntax of inter
inter('clear');
inter('disp');
inter('reset','inter','linearInterMatlab','regularizer','none','theta',1);
inter('disp');
[scheme,parameter] = inter
pause(2);

%% test derivative of schemes
% get all interpolation schemes
files = getFairFiles('interpolation','*Inter*');
for dim=3,
  % initialize dim-dimensional test data
  switch dim,
    case 1,
      % case dim = 1
      TD = [0 0 1 1 0 0]'; m = length(TD);
      omega = [0,m]; 
      xc = getCellCenteredGrid(omega,m);
      mf = 3*m;
      xf = getCellCenteredGrid(omega,mf);
    case 2,
      omega = [0,1,0,2]; m = [13,16];
      xc = getCellCenteredGrid(omega,m);
      yc = reshape(xc,[m,dim]);
      TD = (yc(:,:,1)-0.5).^2 + (yc(:,:,2)-0.75).^2 <= 0.15;
      xf = getCellCenteredGrid(omega,4*m);
      yf = reshape(xf,[4*m,dim]);
    case 3,
      omega = [0,1,0,2,0,1]; m = [13,16,7];
      xc = getCellCenteredGrid(omega,m);
      yc = reshape(xc,[m,dim]);
      TD = (yc(:,:,:,1)-0.5).^2 + (yc(:,:,:,2)-0.75).^2 ...
        + (yc(:,:,:,3)-0.5).^2 <= 0.15;
      xf = getCellCenteredGrid(omega,4*m);
      yf = reshape(xf,[4*m,dim]);
  end;
  for k = 1:length(files),
    fprintf('test [%s]  -- %d of %d\n',files{k},k,length(files))
    inter('reset','inter',files{k});
    TC = inter('coefficients',TD,[],omega,'out',0);
    Tf = inter(TC,omega,xf);
    zc  = xc + 1e-4*randn(size(xc));
    fctn = @(x) inter(TC,omega,x);
    
    figure(2*k-1); clf;
    switch dim, 
      case 1,
        plot(xc,TD,'ro',xf,Tf,'g-');
      case 2,
        plot3(yc(:,:,1),yc(:,:,2),TD,'ro'); hold on;
        sh = surf(yf(:,:,1),yf(:,:,2),reshape(Tf,4*m));
        set(sh,'faceAlpha',0.5);
      case 3,
        subplot(1,2,1);
        imgmontage(TD,omega,m);
        subplot(1,2,2);
        imgmontage(Tf,omega,4*m);
    end;
    title(files{k});
    f = checkDerivative(fctn,zc,'fig',2*k);
    ylabel(files{k});
    pause(2); 
    
  end;
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