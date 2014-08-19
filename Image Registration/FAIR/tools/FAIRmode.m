function FAIRmode(mode)

if ~exist('mode','var'), mode = 'normal'; end;

global PAUSE PRINT DEBUG OUT; 

switch mode,
  case 'normal'
    PAUSE = 1;
    PRINT = 1;
    DEBUG = 0;
    OUT   = 2;
  case 'test'
    PAUSE = 0;
    PRINT = 1;
    DEBUG = 1;
    OUT = 2;
  case 'disp'
    fprintf('[global(PAUSE=%d,PRINT=%d,DEBUG=%d)]\n',PAUSE,PRINT,DEBUG);
    return;
  otherwise,
    error('nyi')
end;
fprintf('[global(PAUSE=%d,PRINT=%d,DEBUG=%d,OUT=%d)]\n',...
  PAUSE,PRINT,DEBUG,OUT);
