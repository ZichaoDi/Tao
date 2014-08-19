% (c) Jan Modersitzki 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function files = getFairFiles(folder,pattern)
% 
% returns all mfiles with patter [pattern] in folder [folder]

function files = getFairFiles(folder,pattern)

FAIRbase = which('FAIRcopyright.m');
if isempty(FAIRbase),
  error('can not locate FAIRcopyright.m')
else
  FAIRbase = fileparts(FAIRbase);
end;

fprintf('FAIRbase:      %s\n',FAIRbase);
fprintf('pattern:       %s\n',pattern);
fprintf('folder+add-on: %s\n',folder);

files  = {};
files1 = lookup(fullfile(FAIRbase,folder),pattern);
if iscell(files1),
  files = files1;
end;
files2 = lookup(fullfile(FAIRbase,'..','addons',folder),pattern);
if iscell(files2),
  files = {files{:},files2{:}};
end;


function files = lookup(folder,pattern)

if ~exist(folder,'dir'),
  files = -1;
  fprintf('folder [%s] does not exists!\n',folder);
  return
end;

addpath(folder)
files = dir(fullfile(folder,[pattern,'.m']));
fprintf('found %d file(s)\n',length(files));
files = {files(:).name};
for k=1:length(files), files{k} = files{k}(1:end-2); end;

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