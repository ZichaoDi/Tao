% (c) Jan Modersitzki 2010/12/23, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% FAIR-BIG-Tutorials: LANDMARK BASED REGISTRATION
%
%
% E5_2D_affine           - use affine     transformation
% E5_2D_quadratic        - use quadratic  transformation
% E5_2D_TPS              - use TPS  transformation
% E5_Hands_TPS           - TPS in short

clc, clear, close all, help(mfilename);

jobs = {
  'E5_2D_affine',    'use affine     transformation',...
  'E5_2D_quadratic', 'use quadratic  transformation',...
  'E5_2D_TPS',       'use TPS  transformation',...
  'E5_Hands_TPS',   'TPS in short'
};

for k=1:length(jobs)/2,
  runThis(jobs{2*k-1},jobs{2*k});
end;

fprintf('\n<%s> done!\n',mfilename); 
