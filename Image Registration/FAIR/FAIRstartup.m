% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% setup toolbox path

function FAIRstartup

fprintf('\n%% %s\n',char('='*ones(78,1)));
fprintf('FAIR: Flexible Algorithms for Image Registration\n');
fprintf('(c) Jan Modersitzki  -- 2011-01-11\n');
fprintf('\n%% %s\n',char('='*ones(78,1)));
fprintf('Set path on [%s], pwd is [%s]\n',computer,pwd);

FAIRpath = fileparts(which('FAIRcopyright.m'));

str = sprintf('function value=FAIRpath; value=''%s'';',FAIRpath);
fprintf('  - %s\n',str);
f = fopen(fullfile('tools','FAIRpath.m'),'w'); fprintf(f,'%s\n',str); fclose(f);

fairTemp = fullfile(FAIRpath,'temp');
if ~exist(fairTemp,'dir'), mkdir(fairTemp);  end;

folder = dir(FAIRpath);
fprintf('  - addpath %s\n',pwd);
addpath(pwd);
for i=1:length(folder)
  if folder(i).isdir ...
      && ~strcmp(folder(i).name(1),'.') ...
      && ~(folder(i).name(1) == '#'),
    f = fullfile(FAIRpath,folder(i).name);
    fprintf('  - addpath %s\n',f);
    addpath(f);
  end;
end;
fprintf('%% %s\n',char('='*ones(78,1)));
folder = dir([FAIRpath filesep 'apps']);
if not(isempty(folder)),
    fprintf('Add FAIR apps to path\n');
    for i=1:length(folder)
        if folder(i).isdir ...
                && ~strcmp(folder(i).name(1),'.') ...
                && ~(folder(i).name(1) == '#'),
            f = fullfile([FAIRpath filesep 'apps'],folder(i).name);
            fprintf('  - addpath %s\n',f);
            addpath(genpath(f));
        end;
    end;
    fprintf('%% %s\n',char('='*ones(78,1)));
end
FAIRmode('normal')
fprintf('\n[new to FAIR? type: help tutorials and/or run BigTutorial2D]\n\n');

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
