% (c) Jan Modersitzki 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Test script files

function testScriptFiles(folder)

if ~exist('folder','var'), folder = 'data';  end;

clc
mfiles = dir(fullfile(folder,'*.m'));
J = find(folder == filesep);
if ~isempty(J),
  folder = folder(J+1:end)
end;
diaryFile = [folder,'.log'];

fprintf('%s of %s, %d mfiles, log -> %s\n',...
  mfilename,folder,length(mfiles),diaryFile)
FAIRpause(2);


if exist(diaryFile,'file'), delete(diaryFile); end;
log = fopen(diaryFile,'w');

for j=1:length(mfiles)
  file = mfiles(j).name;
  p = fopen(file,'r'); contents = fread(p,'char')'; fclose(p);
  if findstr(contents,'function'),
    fprintf('file %d/%d - %s - is propably a function and thus not called\n',...
      j,length(mfiles),file);
    fprintf(log,'file %d/%d - %s - is propably a function and thus not called\n',...
      j,length(mfiles),file);
  else
    fprintf('try %d/%d - %s\n',...
      j,length(mfiles),file);
    
    try
      [p,file] = fileparts(file);
      eval(file);
      fprintf(log,'%-50s OK\n',file)
    catch 
      fprintf(log,'-- %45s not OK\n',file);      
    end;
    
  end;
end;

fclose(log);
return;

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