function FAIRmessage(str,c)

if ~exist('c','var'), c = '=';  end;

fprintf('>> %s  [ %s ]  % s\n',...
  char(ones(1,10)*c),str,char(ones(1,60-length(str))*c));
