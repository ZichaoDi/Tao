% -----------------------------------------------------------------------------------
% (c) Jan Modersitzki 2011-01-11, see FAIR.2 and FAIRcopyright.m.
% Test FAIR package
%
% -----------------------------------------------------------------------------------
% unified framework for diary files
function FAIRdiary(varargin);
caller = dbstack;               % identify the name of the calling function
caller = caller(min(2,length(caller))).name;
File   = @(str) fullfile(FAIRpath,'temp',[date,'-',str,'.asc']);

diaryFile = File(caller);
mode = 'on';
if nargin == 1,
  if strcmp(varargin{1},'off'),
    mode = 'off';
  elseif ~strcmp(varargin{1},'on'),
    diaryFile = File(varargin{1});
  end;
end;

if strcmp(mode,'on'),
  if exist(diaryFile), delete(diaryFile); end;
  builtin('diary',diaryFile);
  
  fprintf('\n%% %s\n',char('='*ones(78,1)));
  fprintf('FAIR: 2011-01-11 create log for %s\n',caller)
  fprintf('      ==> write log to %s\n',...
    diaryFile(findstr(diaryFile,'temp'):end));
  fprintf('%% %s\n',char('='*ones(78,1)));
else
  fprintf('      ==> close log %s\n',diaryFile);
  builtin('diary','off');
  fprintf('%% %s\n',char('='*ones(78,1)));
end;


