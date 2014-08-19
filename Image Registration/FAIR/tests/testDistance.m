% (c) Jan Modersitzki 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Test DISTANCE MEASURE Tools

toolbox = fileparts(which('distance'))
testStart

debit = {
  'contents.m'
  'distance.m'
  'SSD.m'
  'SSDmex.m'
  'SSDmexC.cpp'
  'SSDweighted.m'
  'NCC.m'
  'NGF.m'
  'NGFcross.m'
  'NGFdot.m'
  'NGFmex.m'                                
  'NCCmex.m'
  'NCCmexC.cpp'
  'NGFdotMexC.cpp'                           
  'rhoSpline.m'
  'rhoSplineC.cpp'
  'MI.m'
  'MIcc.m'
  'MIspline.m'
  };

checkToolbox(toolbox,debit,'EDIT',0,'RUN',0,'VERIFY',0,'KEYBOARD',0); 


%% test syntax of inter
distance('clear');
distance('disp');
distance('reset','distance','MIcc',...
  'tol',1e-7,'minT',0,'maxT',256,'nT',60,'minR',0,'maxR',256,'nR',60);
distance('disp');
[scheme,parameter] = distance
pause(2);


%%
fprintf(' ... %s\n','generate hand data and setup transformation model:')
setupHandData
inter('reset','inter','splineInter','regularizer','moments','theta',1e0);
level = 4; omega = MLdata{level}.omega; m = MLdata{level}.m;
[T,R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega(1,:),'out',0);
xc  = getCellCenteredGrid(omega(1,:),m);
Rc  = inter(R,omega(end,:),xc);
FAIRfigure(2);
subplot(1,2,1); viewImage2Dsc(T,omega,m,'title','template','colormap','gray')
subplot(1,2,2); viewImage2Dsc(R,omega,m,'title','reference','colormap','gray')

%% run over several distances
distances = {'SSD','NCC','MIspline','MIcc','NGFdot'};
clear fig; l = 0;
% ------------------- RUN over all distance measures -------------------------------------
for k=1:length(distances);

  distance('reset','distance',distances{k});
  switch distance,
    case {'SSD','NCC'},
    case 'MIcc',
      %setup default parameter
      distance('set','tol',1e-7,...
        'minT',0,'maxT',256,'nT',60,'minR',0,'maxR',256,'nR',40);
      distance('disp');
    case 'NGFdot',
      %setup default parameter
      distance('set','edge',100);
      distance('disp');
  end;


  for type = 1:2,
    fctn = @(x) parseDistance(type,T,Rc,omega(end,:),m,x);
    l = l+1;
    fig(l)=checkDerivative(fctn,xc);
    ylabel([distance,'-',int2str(type)])
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