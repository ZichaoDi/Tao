% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% FAIR-BIG-Tutorials: DISTANCES
% see also
%   E7_Hands_SSDvsRotation   - SSD versus rotation (hands)
%   E7_PETCT_SSDvsRotation   - SSD versus rotation (PET/CT)
%   E7_PETCT_MIvsRotation    - MI versus rotation  (PET/CT)
%   E7_US_MIvsRotation       - MI versus rotation  (US/-US) 
%   E7_basic                 - distances versus rotation 
%   E7_extended              - distances vs rotation (ext)
%   E7_SSDforces             - show SSD forces
%   E6_PETCT_MLPIR           - MLPIR, PETCT, affine, SSD
%   E9_PETCT_MLIR_NGF_mbElas - MLIR, PETCT, NGF, elastic, matrix based
%   E9_HNSP_MLIR_TR          - MLIR, HNSP, elastic, matrix free, TrustRegion


clc, clear, close all, help(mfilename);

jobs = {
%   'E7_Hands_SSDvsRotation',   'SSD versus rotation (hands)',...
%   'E7_PETCT_SSDvsRotation',   'SSD versus rotation (PET/CT)',...  
  'E7_PETCT_MIvsRotation',    'MI versus rotation  (PET/CT)',...   
%   'E7_US_MIvsRotation',       'MI versus rotation  (US/-US) ',...  
%   'E7_basic',                 'distances versus rotation ',...    
%   'E7_extended',              'distances vs rotation (ext)',... 
%   'E7_SSDforces',             'show SSD forces',... 
%   'E6_PETCT_MLPIR',           'MLPIR, PETCT, affine, SSD',...
%   'E9_PETCT_MLIR_NGF_mbElas', 'MLIR, PETCT, NGF, elastic, matrix based',...
%   'E9_HNSP_MLIR_TR',          'MLIR, HNSP, elastic, matrix free, TrustRegion'
};

for k=1:length(jobs)/2,
  runThis(jobs{2*k-1},jobs{2*k});
end;

fprintf('\n<%s> done!\n',mfilename); 
