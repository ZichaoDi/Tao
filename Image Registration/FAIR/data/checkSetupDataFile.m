% (c) Jan Modersitzki 2010/12/25, see FAIR.2 and FAIRcopyright.m.
% Flexible Algorithms for Image Registration, SIAM, 2009
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% If data has already been saved, just load it and initializes parameterization

caller  = dbstack;               % identify the name of the calling function
caller  = caller(min(length(caller),2)).name;
outfile = fullfile(FAIRpath,'temp',[caller,'.mat']);
OK     = exist(outfile,'file');
if ~OK, 
  return; % matfile has to be generated
end;

fprintf('%s: load(%s); and initialize modules\n',caller,outfile(length(FAIRpath)+1:end));
load(outfile);

% initialize toolbox
if ~isempty(who(outfile,'viewOptn')),  viewImage('reset',viewOptn{:});  end;
if ~isempty(who(outfile,'intOptn')),   inter('reset',intOptn{:});       end;
if ~isempty(who(outfile,'traOptn')),   trafo('reset',traOptn{:});       end;
if ~isempty(who(outfile,'disOptn')),   distance('reset',disOptn{:});    end;
if ~isempty(who(outfile,'regOptn')),   regularizer('reset',regOptn{:}); end;
global OUT,
if OUT>1,
  reportStatus
end;

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