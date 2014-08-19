% (c) Jan Modersitzki 2010/12/23, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% FAIR-BIG-Tutorials: TRANSFORMATIONS
% see also
%   E4_US_trafo       - plain rotation of 2D US image
%   E4_US_trafos      - transformations of 2D US image
%   E4_US_rotation    - rotate 2D US image, (+grid)
%   E4_US             - various transformations for 2D US image

clc, clear, close all, help(mfilename);

jobs = {
  'E4_US_trafo',       'plain rotation of 2D US image',...
%   'E4_US_rotation',    'rotate 2D US image, (+grid)',...
%   'E4_US_trafos',      'various transformations for 2D US image'
};

for k=1:length(jobs)/2,
  runThis(jobs{2*k-1},jobs{2*k});
end;

fprintf('\n<%s> done!\n',mfilename); 
