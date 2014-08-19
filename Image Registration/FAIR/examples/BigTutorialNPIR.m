% Tutorial for FAIR: NON-PARAMETRIC IMAGE REGISTRATION
% (c) Jan Modersitzki 2009/04/07, see FAIR.2 and FAIRcopyright.m.
% \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
%
% E9_Hands_NPIRmb_GN:       NPIR, hands, elastic, matrix based
% E9_Hands_NPIRmf_GN:       NPIR, hands, elastic, matrix free
% E9_Hands_NPIRmf_TR_nopre: NPIR, hands, Trust-Region, plain
% E9_Hands_NPIRmf_TR_pcg:   NPIR, hands, Trust-Region, MG-preconditioned
% E9_PETCT_MLIR:            MLIR, PETCT, NGF, elastic, matrix based
% E9_HNSP_MLIR_TR:          MLIR, HNSP, elastic, matrix free, TrustRegion
% E9_Hands_MLIRmb:          MLIR, hands, elastic, matrix based
% E9_Hands_MLIRmf:          MLIR, hands, elastic, matrix free
% E9_Hands_MLIR_mbCurv:     MLIR, hands, curvature, matrix based
% E9_Hands_MLIR_mfCurv:     MLIR, hands, curvature, matrix free
% E9_HNSP_MLIR_SSD_mbElas:  MLIR, HNSP, elastic, matrix based
% E9_HNSP_MLIR_SSD_mfElas:  MLIR, HNSP, elastic, matrix free
% E9_HNSP_MLIR_mbCurv:      MLIR, HNSP, curvature, matrix based
% E9_HNSP_MLIR_mfCurv:      MLIR, HNSP, curvature, matrix free

clc, clear, close all, help(mfilename);

% NPIR: Non-Parametric Image Registration (fixed level)
% runThis('E9_Hands_NPIRmb_GN','NPIR, hands, elastic, matrix based');
% runThis('E9_Hands_NPIRmf_GN','NPIR, hands, elastic, matrix free');
% runThis('E9_Hands_NPIRmf_TR_nopre','NPIR, hands, Trust-Region, plain');
% runThis('E9_Hands_NPIRmf_TR_pcg','NPIR, hands, Trust-Region, MG-preconditioned');
% 
% % MLIR: MultiLevel Non-Parametric Image Registration
% runThis('E9_PETCT_MLIR_NGF_mbElas','MLIR, PETCT, NGF, elastic, matrix based');
% runThis('E9_HNSP_MLIR_TR','MLIR, HNSP, elastic, matrix free, TrustRegion');
% %%
% runThis('E9_Hands_MLIR_SSD_mbElas','MLIR, Hands, elastic, matrix based');
% runThis('E9_Hands_MLIR_SSD_mfElas','MLIR, Hands, elastic, matrix free');
runThis('E9_Hands_MLIR_SSD_mbCurv','MLIR, Hands, curvature, matrix based');
runThis('E9_Hands_MLIR_SSD_mfCurv','MLIR, Hands, curvature, matrix free');

runThis('E9_HNSP_MLIR_SSD_mbElas','MLIR, HNSP, elastic, matrix based');
runThis('E9_HNSP_MLIR_SSD_mfElas','MLIR, HNSP, elastic, matrix free');
runThis('E9_HNSP_MLIR_SSD_mbCurv','MLIR, HNSP, curvature, matrix based');
runThis('E9_HNSP_MLIR_SSD_mfCurv','MLIR, HNSP, curvature, matrix free');

fprintf('\n<%s> done!\n',mfilename); 
