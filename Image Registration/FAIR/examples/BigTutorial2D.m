% Tutorial for FAIR: 2D EXAMPLES
% (c) Jan Modersitzki 2009/04/07, see FAIR.2 and FAIRcopyright.m.
% \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
%
% GENERAL: 
% E2_setupHandData:        initialize 2D data for use in FAIR
% E2_setupBrainData:       initialize 3D data for use in FAIR
% E2_viewImage:            usinig the image viewer
% E3_Hands_ij2xy:          transform left handed to right handed data
% E3_US_getMultilevel:     generate MultiLevel representation of data
% E4_US_rotation:          rotate an US image wrt. domain center
% E4_US_trafos:            same as above, but various transformations
%
% PARAMETRIC IMAGE REGISTRATION
% E5_Hands_TPS:            Thin-Plate-Spline Registration
% E6_Hands_PIR_SD:         PIR, hands, rotations, SSD, Steepest Descent
% E6_Hands_PIR_GN:         PIR, hands, rotations, SSD, Gauss-Newton
%
% MULTILEVEL PARAMETRIC IMAGE REGISTRATION
% E6_Hands_MLPIR:          MLPIR, hands, affine, SSD
% E6_Hands_MLPIR_pause:    like above, including pauses
% E6_HNSP_MLPIR_reg:       MLPIR, HNSP, splines, regularizes
% E6_PETCT_MLPIR:          MLPIR, PETCT, affine, SSD
%
% NON-PARAMETRIC IMAGE REGISTRATION
% E9_Hands_NPIRmb_GN:       NPIR, hands, elastic, matrix based
% E9_Hands_NPIRmf_GN:       NPIR, hands, elastic, matrix free
% E9_Hands_NPIRmf_TR_nopre: NPIR, hands, Trust-Region, plain
% E9_Hands_NPIRmf_TR_pcg:   NPIR, hands, Trust-Region, MG-preconditioned
%
% MULTILEVEL NON-PARAMETRIC IMAGE REGISTRATION
% E9_PETCT_MLIR_NGF_mbElas: MLIR, PETCT, NGF, elastic, matrix based
% E9_HNSP_MLIR_TR:          MLIR, HNSP, elastic, matrix free, TrustRegion
% E9_Hands_MLIR_SSD_mbElas: MLIR, hands, elastic, matrix based
% E9_Hands_MLIR_SSD_mfElas: MLIR, hands, elastic, matrix free
% E9_Hands_MLIR_SSD_mbCurv: MLIR, hands, curvature, matrix based
% E9_Hands_MLIR_SSD_mfCurv: MLIR, hands, curvature, matrix free
% E9_HNSP_MLIR_SSD_mbElas:  MLIR, HNSP, elastic, matrix based
% E9_HNSP_MLIR_SSD_mfElas:  MLIR, HNSP, elastic, matrix free
% E9_HNSP_MLIR_SSD_mbCurv:  MLIR, HNSP, curvature, matrix based
% E9_HNSP_MLIR_SSD_mfCurv:  MLIR, HNSP, curvature, matrix free

clc, clear, close all, help(mfilename);

% general stuff
runThis('E2_setupHandData','initialize 2D data for use in FAIR');
runThis('E2_setupBrainData','initialize 3D data for use in FAIR');
runThis('E2_viewImage','usinig the image viewer');
runThis('E3_Hands_ij2xy','transform left handed to right handed data');
runThis('E3_US_getMultilevel','generate MultiLevel representation of data');
runThis('E4_US_rotation','rotate an US image wrt. domain center');
runThis('E4_US_trafos','same as above, but various transformations');
runThis('E5_Hands_TPS','Thin-Plate-Spline Registration');

% Parametric Image Registration  (fixed level)
runThis('E6_Hands_PIR_SD','PIR, hands, rotations, SSD, Steepest Descent');
runThis('E6_Hands_PIR_GN','PIR, hands, rotations, SSD, Gauss-Newton');

% MultiLevel Parametric Image Registration
runThis('E6_Hands_MLPIR','MLPIR, hands, affine, SSD');
runThis('E6_Hands_MLPIR_pause','like above, including pauses');
runThis('E6_HNSP_MLPIR_reg','MLPIR, HNSP, splines, regularizes');
runThis('E6_PETCT_MLPIR','MLPIR, PETCT, affine, SSD');

% Non-Parametric Image Registration (fixed level)
runThis('E9_Hands_NPIRmb_GN','NPIR, hands, elastic, matrix based');
runThis('E9_Hands_NPIRmf_GN','NPIR, hands, elastic, matrix free');
runThis('E9_Hands_NPIRmf_TR_nopre','NPIR, hands, Trust-Region, plain');
runThis('E9_Hands_NPIRmf_TR_pcg','NPIR, hands, Trust-Region, MG-preconditioned');

% MultiLevel Non-Parametric Image Registration
runThis('E9_PETCT_MLIR_NGF_mbElas','MLIR, PETCT, NGF, elastic, matrix based');
runThis('E9_HNSP_MLIR_TR','MLIR, HNSP, elastic, matrix free, TrustRegion');

runThis('E9_Hands_MLIR_SSD_mbElas','MLIR, hands, elastic, matrix based');
runThis('E9_Hands_MLIR_SSD_mfElas','MLIR, hands, elastic, matrix free');
runThis('E9_Hands_MLIR_SSD_mbCurv','MLIR, hands, curvature, matrix based');
runThis('E9_Hands_MLIR_SSD_mfCurv','MLIR, hands, curvature, matrix free');

runThis('E9_HNSP_MLIR_SSD_mbElas','MLIR, HNSP, elastic, matrix based');
runThis('E9_HNSP_MLIR_SSD_mfElas','MLIR, HNSP, elastic, matrix free');
runThis('E9_HNSP_MLIR_SSD_mbCurv','MLIR, HNSP, curvature, matrix based');
runThis('E9_HNSP_MLIR_SSD_mfCurv','MLIR, HNSP, curvature, matrix free');

fprintf('\n<%s> done!\n',mfilename); 
