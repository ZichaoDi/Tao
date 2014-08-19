% wrapper for tests while calling a tutorial
function run = runThis(file,str,varargin)

%global PAUSE
PAUSE=1;
% return
test = 1;
for k=1:2:length(varargin), % overwrites default parameter 
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

fprintf('%%----- %-24s %-30s',[file,':'],str);

if PAUSE,
  run = input('<0,1> : ');
else
  run = test;
  fprintf('\n\n global(PAUSE==0) -> run(TEST=%d)\n',test)
  clear = @()         fprintf('clear has been disabled for auto-testing\n')
  clc   = @()         fprintf('clc   has been disabled for auto-testing\n')
  pause = @(varargin) fprintf('pause has been disabled for auto-testing\n')

end;

if ~run,            return; end;
if isempty(file),   return; end;
eval(file);
