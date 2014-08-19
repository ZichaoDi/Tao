% Tutorial for FAIR: 3D EXAMPLES
% (c) Jan Modersitzki 2009/04/07, see FAIR.2 and FAIRcopyright.m.
% \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
%
% E2_viewImage3D:          initialize and visualize 3D data
% E2_setupBrainData:       initialize 3D data for use in FAIR
% E9_3Dbrain_GN:           MLIR for brain data (Gauss-Newton)
% E9_3Dbrain_TR:           MLIR for brain data (Trust-Region)
% E9_3Dbrain_lBFGS:        MLIR for brain data (lBFGS)

clc, clear, close all, help(mfilename);

runThis('E2_viewImage3D','initialize and visualize 3D data');
runThis('E2_setupBrainData','initialize 3D data for use in FAIR');
runThis('E9_3Dbrain_GN','MLIR for brain data (Gauss-Newton)');
runThis('E9_3Dbrain_TR','MLIR for brain data (Trust-Region)');
runThis('E9_3Dbrain_lBFGS','MLIR for brain data (lBFGS)');

fprintf('\n<%s> done!\n',mfilename); 
