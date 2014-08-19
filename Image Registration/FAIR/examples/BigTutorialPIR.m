% (c) Jan Modersitzki 2010/12/23, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% FAIR-BIG-Tutorials: PARAMETRIC IMAGE REGISTRATION
% see also
%   E7_Hands_SSDvsRotation  - SSD versus rotation (hands)
%%
%   E6_Hands_PIR_SD         - PIR, hands, rotations, SSD, Steepest Descent
%   E6_Hands_PIR_GN         - PIR, hands, rotations, SSD, Gauss-Newton
%%
%   E6_Hands_MLPIR_pause    - MLPIR, hands, affine, SSD, including pauses
%   E6_Hands_MLPIR          - MLPIR, hands, affine, SSD
%   E6_HNSP_MLPIR_reg       - MLPIR, HNSP, splines, regularizes
%   E6_PETCT_MLPIR          - MLPIR, PETCT, affine, SSD

clc, clear, close all, help(mfilename);

jobs = {
%    'E7_Hands_SSDvsRotation', 'SSD versus rotation (hands)',...
   'E7_Lena_SSDvsRigid', 'SSD versus rigid (Lena)',...
%   'E6_Hands_PIR_SD','PIR, hands, rotations, SSD, Steepest Descent',...
%   'E6_Hands_PIR_GN','PIR, hands, rotations, SSD, Gauss-Newton',...
%   'E6_Hands_MLPIR_pause','MLPIR, hands, affine, SSD, including pauses',...
%   'E6_Hands_MLPIR','MLPIR, hands, affine, SSD',...
%   'E6_HNSP_MLPIR_reg','MLPIR, HNSP, splines, regularizes',...
%   'E6_PETCT_MLPIR','MLPIR, PETCT, affine, SSD'
};

for k=1:length(jobs)/2,
  runThis(jobs{2*k-1},jobs{2*k});
end;

fprintf('\n<%s> done!\n',mfilename); 
 
