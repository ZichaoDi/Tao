% (c) Jan Modersitzki 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% test the setup of various examples

toolbox = fileparts(which('setupHandData'));
testStart


debit = {
  'get2Ddata.m'
  'setupUSData.m'
  'setup2Ddisc2CData.m'
  'setup2DEPIData.m'
  'setup2DGaussianData.m'
  'setup2DGaussianNoisyData.m'
  'setup2DMPsquare.m'
  'setupHNSPData.m'
  'setupHandData.m'
  'setupMRIData.m'
  'setupMyData.m'
  'setupPETCTData.m'
  'setup3DboxData.m'
  'setup3DbrainData.m'
  'setup3DkneeData.m'
  'setup3DmiceData.m'
  'setup3DphantomData.m'
  'contents.m'
  'checkSetupDataFile.m'
  'c.jpg'
  'circle.jpg'
  'EPIslice-R.jpg'
  'EPIslice-T.jpg'
  'Gauss-R.jpg'
  'Gauss-T.jpg'
  'Gauss-R-noisy.jpg'
  'Gauss-T-noisy.jpg'
  'hands-R.jpg'
  'hands-T.jpg'
  'HNSP-R.jpg'
  'HNSP-T.jpg'
  'MPSquare-R.jpg'
  'MPSquare-T.jpg'
  'MRIhead-R.jpg'
  'MRIhead-T.jpg'
  'PET-CT-CT.jpg'
  'PET-CT-PET.jpg'
  'US Vibe Heldmann.jpg'
  'US.jpg'
  'brain2D.mat'
  'brain3D.mat'
  'createPhilipsKnees.mat'
  'mice3D.mat'
  'phantom3D.mat'
  };

checkToolbox(toolbox,debit,'EDIT',0,'RUN',0);

% list of mandatory variables
list = {'dataT','omega','MLdata','m','viewOptn'};

mFiles = cell(debit);
j = find(strcmp(mFiles,[mfilename,'.m']));
mFiles(j) = []; % do not run the calling file recursively!

% RUN FILE
for k=1:length(mFiles),
  
  [p,f,e] = fileparts(mFiles{k}(1:end));
  
  if strcmp(e,'.m'),
    fprintf('\n\n\n\n\n  - run  %-3d of %d, %20s\n\n',...
      k,length(mFiles),mFiles{k});
    mat = which([mFiles{k}(1:end-2),'.mat']);
    if ~isempty(mat), delete(mat); end;
    run(mFiles{k}(1:end-2));
    if ~isempty(findstr(mFiles{k},'setup')) && ~inMatFile(mFiles{k}(1:end-2),list),
      error('missing variables');
    end;
    pause(2)
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