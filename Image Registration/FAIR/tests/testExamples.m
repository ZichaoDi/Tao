% (c) Jan Modersitzki 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Test Examples and such schemes

toolbox = fileparts(which('BigTutorial2D'));

testStart

debit = {
  %{
  % {
  'contents.m'
  'E2_setupBrainData.m'
  'E2_setupHandData.m'
  'E2_viewImage.m'
  'E2_viewImage3D.m'
  % }
  % {
  'E3_1D_basics.m'
  'E3_1D_derivatives.m'
  'E3_1D_scale.m'
  'E3_2D_basics.m'
  'E3_2D_derivative.m'
  'E3_2D_generic.m'
  'E3_2D_scale.m'
   % }
  % {
  'E3_getCellCenteredGrid.m'
  'E3_ij2xy3d.m'
  'E3_Hands_ij2xy.m'
  'E3_multilevel.m'
  'E3_US_getMultilevel.m'
  'E3_matlabInterpolation1D.m'
  'E3_linearInterpolation1D.m'
  'E3_linearInterpolation2D.m'
  'E3_bsplines.m'
  'E3_checkDerivative.m'
  'E3_splineInterpolation1D.m'
  'E3_splineInterpolation2D.m'
  'E3_truncatedSplineInterpolation1D.m'
  'E3_MS_splineInterpolation1D.m'
  'E3_MS_splineInterpolation2D.m'
  'E3_MS_splineInterpolation2Dext.m'
  'E3_interpolation2D.m'
  % }
  % {
  'E4_Affine2D.m'
  'E4_Affine2Dplain.m'
  'E4_Bizarr.m'
  'E4_Rigid2D.m'
  'E4_Rigid2Dplain.m'
  'E4_SplineTransformation2D.m'
  'E4_Translation2D.m'
  'E4_Translation2Dplain.m'
  'E4_US_rotation.m'
  'E4_US_trafo.m'
  'E4_US_trafos.m'
  % }
  % {
  'E5_2D_affine.m'
  'E5_2D_quadratic.m'
  'E5_Hands_TPS.m'
  'E5_2D_TPS.m'
  'E5_linear.m'
  'E5_quadratic.m'
  'E5_TPS.m'
  % }
  % {
  'E6_quadrature_Spline1D.m'
  'E6_quadrature_Spline2D.m'
  'E6_quadrature_Gaussian2D.m'
  'E6_quadrature_SSD2D.m'
  %%   
  'E6_HNSP_SSD_rotation2D_level4.m'
  'E6_HNSP_SSD_rotation2D_level8.m'
  'E6_HNSP_SSD_translation2D_level4.m'
  'E6_HNSP_SSD_translation2D_level8.m'
  'E6_HNSP_SSD_translation2D_level4_spline.m'
  'E6_HNSP_SSD_translation2D_level8_spline.m'
  % }  
   % {
  'E6_FAIRplots.m'
  'E6_HNSP_PIR_SSD_rotation2D_level4.m'
  'E6_HNSP_PIR_SSD_rotation2D_level7.m'
  'E6_HNSP_PIR_SSD_rigid2D_level4.m'
  'E6_HNSP_PIR_SSD_rigid2D_level7.m'
  'E6_HNSP_PIR_SSD_affine2D_level5.m'
  'E6_HNSP_PIR_SSD_spline2D_level5.m'
  %%
  'E6_HNSP_RPIR.m'
  'E6_HNSP_PIR_scale.m'
  'E6_HNSP_PIR_NelderMead.m'
  'E6_HNSP_PIR_GN.m'
  'E6_Hands_PIR_GN.m'
  'E6_Hands_PIR_SD.m'
  %%
  'E6_HNSP_MLPIR_SSD_rotation2D.m'
  'E6_HNSP_MLPIR_SSD_rigid2D.m'
  'E6_HNSP_MLPIR_SSD_affine2D.m'
  'E6_HNSP_MLPIR_reg.m'
  'E6_Hands_MLPIR.m'
  'E6_Hands_MLPIR_pause.m'
  'E6_PETCT_MLPIR.m'
  % }
  % {
  'E7_histogram1D.m'
  'E7_histogram1D_ext.m'
  'E7_SSDforces.m'
  'E7_HNSP_SSD_forces.m'
  'E7_Hands_SSDvsRotation.m'
  'E7_Hands_distance_rotation.m'
  'E7_Hands_distance_rotation_ext.m'
  'E7_PETCT_SSDvsRotation.m'
  'E7_PETCT_MIvsRotation.m'
  'E7_US_MIvsRotation.m'
  'E7_PETCT_MLPIR.m'
  'E7_PETCT_MLPIR_ext.m'
  'E7.m'
  'E7_basic.m'
  'E7_extended.m'
  % }
  % {
  'E8_Belastic.m'
  'E8_checkOperations.m'
  'E8_forces.m'
  'E8_forces_curvature.m'
  'E8_matrices.m'
  'E8_matrices_curvature.m'
  'E8_regularizationMB.m'
  'E8_regularizationMF.m'
  'E8_hyperElasticRegularizationMF.m'
  'E8_hyperElasticRegularizationMB.m'
  'E8_setup.m'  
  'E8_setup_curvature.m'
  %}
  % {
  'E9_HNSP_NPIR.m'
  'E9_HNSP_NPIR_pre.m'
  'E9_Hands_NPIR.m'
  'E9_Hands_NPIR_MI_mbElas_BFGS.m'
  'E9_Hands_NPIR_pre.m'
  'E9_Hands_NPIRmb_GN.m'
  'E9_Hands_NPIRmf_GN.m'
  'E9_Hands_NPIRmf_TR_nopre.m'  %%% TODO: Trust-Region
  'E9_Hands_NPIRmf_TR_pcg.m'    %%% TODO: Trust-Region
  'E9_Hands_affine.m'
  %}
  % {
  'E9_Hands_MLIR_SSD_mbElas.m'
  'E9_Hands_MLIR_SSD_mfElas.m'
  'E9_Hands_MLIR_SSD_mbCurv.m'
  'E9_Hands_MLIR_SSD_mfCurv.m'
  'E9_Hands_MSIR.m'
  %}
  % {
  'E9_HNSP_MLIR_SSD_mbElas.m'
  'E9_HNSP_MLIR_SSD_mfElas.m'
  'E9_HNSP_MLIR_SSD_mbCurv.m'
  'E9_HNSP_MLIR_SSD_mfCurv.m'
  'E9_HNSP_MLIR_TR.m'
  %}
  % {
  'E9_MRIhead_MLIR_SSD_mbElas.m'
  'E9_MRIhead_MLIR_MI_mbElas.m'
  'E9_MRIhead_MLIR_NGF_mbElas.m'
  'E9_MRIhead_MLIRlBFGS_MI_mfElas.m' 
  'E9_MRIhead_MLIRlBFGS_NGF_mbElas.m'
  'E9_MRIhead_MLIRlBFGS_NGF_elas.m'
   %}
  % {
  'E9_PETCT_MLIR_NGF_mbElas.m'
  'E9_PETCT_MLIR_NGF_mfElas.m'
  'E9_PETCT_MLIR_NGF_mbCurv.m'
  %
  'E9_PETCT_MLIRlBFGS_MI_mbCurv.m'
  'E9_PETCT_MLIRlBFGS_MI_mbElas.m'
  'E9_PETCT_MLIRlBFGS_NGF_mbCurv.m'
  'E9_PETCT_MLIRlBFGS_NGF_mbElas.m'
  %
  'E9_3Dbrain_GN.m'
  'E9_3Dbrain_TR.m'
  'E9_3Dbrain_lBFGS.m'
  %}
  'tutWebsite.m'
  };

checkToolbox(toolbox,debit,'EDIT',0,'VERIFY',0,'KEYBOARD',0,'RUN',1)

testEnd


%{ 
	=======================================================================================
	FAIR: Flexible Algorithms for Image Registration, Version 2011
	Copyright (c): Jan Modersitzki
	Maria-Goeppert-Str. 1a, D-23562 Luebeck, Germany
	Email: jan.modersitzki@mic.uni-luebeck.de
	URL:   http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
	=======================================================================================
	No part of this code may be reproduced, stored in a retrieval system,
	translated, transcribed, transmitted, or distributed in any form
	or by any means, means, manual, electric, electronic, electro-magnetic,
	mechanical, chemical, optical, photocopying, recording, or otherwise,
	without the prior explicit written permission of the authors or their
	designated proxies. In no event shall the above copyright notice be
	removed or altered in any way.

	This code is provided "as is", without any warranty of any kind, either
	expressed or implied, including but not limited to, any implied warranty
	of merchantibility or fitness for any purpose. In no event will any party
	who distributed the code be liable for damages or for any claim(s) by
	any other party, including but not limited to, any lost profits, lost
	monies, lost data or data rendered inaccurate, losses sustained by
	third parties, or any other special, incidental or consequential damages
	arrising out of the use or inability to use the program, even if the
	possibility of such damages has been advised against. The entire risk
	as to the quality, the performace, and the fitness of the program for any
	particular purpose lies with the party using the code.
	=======================================================================================
	Any use of this code constitutes acceptance of the terms of the above statements
	=======================================================================================
%}