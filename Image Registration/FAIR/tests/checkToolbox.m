% (c) Jan Modersitzki 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function checkToolbox(toolbox,debit,varargin)

function checkToolbox(toolbox,debit,varargin)

EDIT     = 1;
RUN      = 1;
MEX      = 1;
VERIFY   = 0;
KEYBOARD = 0;
range    = [];
clear = @()         fprintf('clear has been disabled for auto-testing\n')
clc   = @()         fprintf('clc   has been disabled for auto-testing\n')
pause = @(varargin) fprintf('pause has been disabled for auto-testing\n')

for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if exist(toolbox) & ~ isempty(toolbox),
  cd(toolbox);
end;

mFiles = cell(debit);
credit = dir('.');

fileCheck = 1;

fprintf('\nFILES IN TOOLBOX [%s]\n',toolbox);
for k=1:length(debit),
  j = find(strcmp({credit(:).name},debit{k}));
  if isempty(j),
    fprintf('  -> %-40s NOT in %s\n',debit{k},toolbox);
    fileCheck = 0;
  else
    fprintf('  -> %-40s exists\n',debit{k});
    credit(j) = [];
  end;
end;

if ~fileCheck, error('TOOLBOX INCOMPLETE --- FILES ARE MISSING');  end;

% remove "*.m~"
names = {credit(:).name};
J = [find(cellfun(@isempty,strfind(names,'.m~')) == 0)];
credit(J) = [];
if length(credit) > 0,
  fprintf('\nADDITIONAL FILES IN [%s]\n',toolbox);
  for k=1:length(credit),
    if credit(k).name(1) ~= '.',
      fprintf('  -> %-40s\n',credit(k).name);
    end;
  end;
end;


if MEX,
  fprintf('\nREBUILD MEX-FILES IN TOOLBOX [%s]\n',toolbox);
  cFiles = [dir('*.c');dir('*.C');dir('*.cpp');dir('*.CPP')];
  for k=1:length(cFiles),
    FAIRmake(cFiles(k).name)
  end;
end;


% VERIFY against old version
if VERIFY,
%   oldFolder = fullfile('/Users/Jan/JM/FAIR/Versions/SIAM-FAIR',...
%     toolbox(find(toolbox==filesep,1,'last'):end));
  oldFolder = fullfile('/Users/Jan/JM/FAIR/Versions','SIAM-FAIR');
  fprintf('\nCOMPARE TO PREVIOUS VERSION:\n');
  fprintf('  [%s]\n',oldFolder);
  
  for k=1:length(mFiles),
    fprintf('  - %-3d of %d, %-40s - ',...
      k,length(mFiles),mFiles{k});
    
    % locate oldFile
    checkFile = evalc(sprintf('!find "%s" -name "%s.checked"',oldFolder,mFiles{k}));
    oldFile   = evalc(sprintf('!find "%s" -name "%s"',oldFolder,mFiles{k}));
    if ~isempty(checkFile),
      fprintf(' already marked checked\n');
    elseif ~isempty(oldFile),
      oldFile(end) = [];
      str = sprintf('!opendiff "%s" "%s" ',mFiles{k},oldFile);
      disp(str);
      eval(str);
      keyboard
      str = sprintf('!mv "%s" "%s.checked" ',oldFile,oldFile);
      disp(str);
      eval(str);
    else
      fprintf(' NEW file\n')
      edit('README_new.m');
    end;
  end;
end;



if isempty(range), 
  range = 1:length(debit); 
else
  range(range>length(debit)) = [];
end;
fileList = debit(range)

% EDIT FILE
if EDIT,
  fprintf('\nEDIT M-FILES IN TOOLBOX [%s]\n',toolbox);
  for k=1:length(fileList),
    fprintf('  - edit %-3d of %d, %20s\n',...
      k,length(fileList),fileList{k});
    
    edit(fileList{k});
    
    %if k>10, return; end;
  end;
end;

% pick only .m files
caller = dbstack;               % identify the name of the calling function
j = find(strcmp(fileList,[caller(2).name,'.m']));
J = [find(cellfun(@isempty,strfind(fileList,'.m')));j];
fileList(J) = []; % do not run the calling file recursively!


% RUN FILE
if RUN,
  for k=1:length(fileList),
    fprintf('\n\n\n\n\n  - run  %-3d of %d, %20s\n\n',...
      k,length(fileList),fileList{k}(1:end-2));
    
    run(fileList{k}(1:end-2))
    if KEYBOARD,
      keyboard
    else
      pause(2)
    end;
  end;
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