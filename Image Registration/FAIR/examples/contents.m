% (c) Jan Modersitzki 2010-08-11, see FAIR.2 and FAIRcopyright.m.
% \url{http://www.mic.uni-luebeck.de/people/jan-modersitzki/software/fair.html}
%
% FAIR Examples Toolbox, contents.m: this file
% 
% The examples are grouped in different categories:
% 
%  B     BigTutorial*.m: Tutorials covering a topic, extended interaction
%  E     E?_*.m:         Examples related to Chapter ? in the FAIR book         
% 
%  E1    has no mfiles
% 
%  E2    BASIC CONCEPTS: DATA LOADING AND VISUALIZATION
% 
%  E3    INTERPOLATION
%  3.1   GRIDS AND COORDINATE SYSTEMS
%  3.2   1D INTERPOLATION
%  3.3   2D INTERPOLATION
%  3.4   MULTISCALE INTERPOLATION
%  3.5   MULTILEVEL REPRESENTATION OF DATA
%  3.6   DERIVATIVES
%  3.7   MISC
%
%  E4    TRANSFORMATIONS
% 
%  E5    LANDMARK REGISTRATION
% 
%  E6    PARAMETRIC IMAGE REGISTRATION
%  6.1   QUADRATURE
%  6.2   EXPLORING THE OBJECTIVE FUNCTION (USING SSD)
%  6.3   PIR: PARAMETRIC IMAGE REGISTRATION
%  6.4   RPIR: REGULARIZED PARAMETRIC IMAGE REGISTRATION
%  6.5   MLPIR: MULTILEVEL PARAMETRIC IMAGE REGISTRATION
%  6.6   VISUALIZATION
% 
%  E7    DISTANCE MEASURES
%  7.1   DISTANCE VERSUS PARAMETERS
%  7.2   FORCES
%  7.3   HISTOGRAMS
%  7.4   DISTANCES in MLPIR
% 
%  E8    REGULARIZATION
% 
%  E9    NON-PARAMETRIC IMAGE REGISTRATION
%  9.1   FIXED LEVEL
%  9.2   MULTILEVEL
%  9.3   MULTISCALE
%  9.4   3D DATA
% 
%  U     AS SHOWN IN THE FAIR BOOK
% 
%  S     SETUP DATA
% 
%  X     MISCELLANEOUS           
% ==============================================================================
% 
%  B     BIGTUTORIAL*: TUTORIALS COVERING A TOPIC, EXTENDED INTERACTION
% 
%      - BigTutorial2D               2D examples
%      - BigTutorial3D               3D examples
%      - BigTutorialDistance         aspects of distance measures
%      - BigTutorialInter            interpolation issues
%      - BigTutorialLandmarks        landmark registration
%      - BigTutorialNPIR             non-parametric examples
%      - BigTutorialPIR  	           parametric examples
%      - BigTutorialRegularizer      regularization
%      - BigTutorialTrafo            transformation 
%      - runThis                     generic alert, use to handle the test suite
% 
% ==============================================================================
%      
%  E2    BASIC CONCEPTS: DATA LOADING AND VISUALIZATION
% 
%      - E2_setupBrainData      example for 3D data  ???
%      - E2_setupHandData       example for 2D data   §
%      - E2_viewImage           visualization in 2D  §
%      - E2_viewImage3D         visualization in 2D  §
% 
% ==============================================================================
%      
%  E3    INTERPOLATION
% 
%  3.1   GRIDS AND COORDINATE SYSTEMS 
% 
%      - E3_Hands_ij2xy                 changing to a right handed coordinate system 
%      - E3_getCellCenteredGrid             creates cell-centered grids 
%      - E3_ij2xy3d                     changing to a right handed coordinate system
% 
%  3.2   1D INTERPOLATION 
% 
%      - E3_1D_basics                   basic concepts 
%      - E3_1D_derivatives              derivative example 
%      - E3_1D_scale                    scale example 
%      - E3_bsplines introducing        B-splines
%      - E3_matlabInterpolation1D       explaining usage of MATLAB''s interpolation 
%      - E3_linearInterpolation1D       explaining usage of linear interpolation 
%      - E3_splineInterpolation1D       explaining usage of spline interpolation 
%      - E3_truncatedSplineInterpolation1D   truncated frequency representation 
%      - spline1D                       1D spline implementation (for plots only)
% 
%  3.3   2D INTERPOLATION 
% 
%      - E3_2D_basics                   basic concepts 
%      - E3_2D_derivative               derivative example 
%      - E3_2D_generic                  interpolation, how to use it 
%      - E3_2D_scale                    scale example 
%      - E3_linearInterpolation2D       explaining usage of linear interpolation 
%      - E3_splineInterpolation2D       explaining usage of spline interpolation 
%      - E3_interpolation2D             how to use FAIR interpolation 
%      - splineInterpolation2D          explaining 2D spline interpolation   
% 
%  3.4   MULTISCALE INTERPOLATION 
% 
%      - E3_MS_splineInterpolation1D 
%      - E3_MS_splineInterpolation2D 
% 
%  3.5   MULTILEVEL REPRESENTATION OF DATA 
% 
%      - E3_US_getMultilevel            create multilevel representation of data 
%      - E3_multilevel                  how to use the multi-level representation 
% 
%  3.6   DERIVATIVES
% 
%      - E3_checkDerivative             how to check the derivative of an
%                                       interpolation scheme 
%  3.7   MISCALENEOUS
%
%      - E3_bsplines                    plots a B-spline 
% 
% ==============================================================================
% 
%  E4    TRANSFORMATIONS
% 
%      - E4_Affine2D                affine linear transformations using FAIR trafo
%      - E4_Affine2Dplain           plain affine linear transformations 
%      - E4_Bizarr                  strange transformations
%      - E4_Rigid2D                 rigid transformations using FAIR trafo
%      - E4_Rigid2Dplain            plain rigid transformations  
%      - E4_SplineTransformation2D  spline transformations using FAIR trafo 
%      - E4_Translation2D           translation using FAIR trafo 
%      - E4_Translation2Dplain      plain translation 
%      - E4_US                      rotating an ultrasound image 
%      - E4_US_rotation             modification of the above (skim version) 
%      - E4_US_trafo                more transformations of the ultrasound image 
%      - E4_US_trafos               modification of the above 
% 
% ==============================================================================
% 
%  E5    LANDMARK REGISTRATION
% 
%      - E5_2D_TPS                 thin-plate-spline landmark registration
%      - E5_2D_affine              affine linear landmark registration
%      - E5_2D_quadratic           quadratic landmark registration
%      - E5_Hands_TPS              thin-plate-spline landmark registration 
%      - E5_TPS                    thin-plate-spline landmark registration
%      - E5_linear                 affine linear landmark registration
%      - E5_quadratic              quadratic landmark registration
%      - P5_LM                     generic plotting driver     
% 
% ==============================================================================
% 
%  E6    PARAMETRIC IMAGE REGISTRATION
% 
%  6.1   QUADRATURE
% 
%      - E6_quadrature_Gaussian2D        numerical integration of a 2D Gaussian
%      - E6_quadrature_SSD2D             numerical integration of the SSD distance measure
%      - E6_quadrature_Spline1D          numerical integration of a 1D Spline
%      - E6_quadrature_Spline2D          numerical integration of a 2D Spline
% 
%  6.2   EXPLORING THE OBJECTIVE FUNCTION (USING SSD)
% 
%      - E6_HNSP_SSD_rotation2D_level4            rotations on level=4
%      - E6_HNSP_SSD_rotation2D_level8            rotations on level=8
%      - E6_HNSP_SSD_translation2D_level4         translations on level=4, linear
%      - E6_HNSP_SSD_translation2D_level4_spline  translations on level=4, spline
%      - E6_HNSP_SSD_translation2D_level8         translations on level=8, linear
%      - E6_HNSP_SSD_translation2D_level8_spline  translations on level=8, spline 
% 
%  6.3   PIR: PARAMETRIC IMAGE REGISTRATION
% 
%      - E6_Hands_PIR_GN                          Hands with Gauss-Newton
%      - E6_Hands_PIR_SD                          Hands with Steepest Descent 
%      - E6_HNSP_PIR_GN                           HNSP with Gauss-Newton
%      - E6_HNSP_PIR_NelderMead                   HNSP with Nelder-Mead 
%      - E6_HNSP_PIR_SSD_affine2D_level5          HNSP using affine transformation, level=5
%      - E6_HNSP_PIR_SSD_rigid2D_level4           HNSP using rigid  transformation, level=4
%      - E6_HNSP_PIR_SSD_rigid2D_level7           HNSP using rigid  transformation, level=7
%      - E6_HNSP_PIR_SSD_rotation2D_level4        HNSP using rotation,              level=4
%      - E6_HNSP_PIR_SSD_rotation2D_level7        HNSP using rotation,              level=7
%      - E6_HNSP_PIR_scale                        HNSP using a scale space
%      - E6_HNSP_PIR_spline2D_level5              HNSP using spline transformation, level=5
% 
%  6.4   RPIR: REGULARIZED PARAMETRIC IMAGE REGISTRATION
% 
%      - E6_HNSP_RPIR                             HNSP using a regularized spline transformation
% 
%  6.5   MLPIR: MULTILEVEL PARAMETRIC IMAGE REGISTRATION
% 
%      - E6_Hands_MLPIR                           Hands
%      - E6_Hands_MLPIR_pause                     as above including pauses 
%      - E6_HNSP_MLPIR_SSD_affine2D               HNSP using affine linear transformation
%      - E6_HNSP_MLPIR_SSD_rigid2D                HNSP using rigid transformation
%      - E6_HNSP_MLPIR_SSD_rotation2D             HNSP using rotations
%      - E6_HNSP_MLPIR_reg                        HNSP using regularized splines
%      - E6_PETCT_MLPIR                           PET-CT data
% 
%  6.6   VISUALIZATION
% 
%      - E6_FAIRplots                             how to use FAIRplots
% 
% ==============================================================================
% 
%  E7    DISTANCE MEASURES
% 
%  7.1   DISTANCE VERSUS PARAMETERS
% 
%      - E7_Hands_SSDvsRotation             SSD versus rotations
%      - E7_Hands_distance_rotation         various distances versus rotations
%      - E7_Hands_distance_rotation_ext     extended version of the above
%      - E7_PETCT_MIvsRotation              MI versus rotations
%      - E7_PETCT_SSDvsRotation             SSD versus rotations
%      - E7_US_MIvsRotation                 MI versus rotations
%      - E7_basic                           various distances versus rotations
%      - E7_extended                        extended version of the above
%      - E7                                 various distances versus rotations and translations
% 
%  7.2   FORCES
% 
%      - E7_SSDforces                       displays SSD force field
%      - E7_HNSP_SSD_forces                 displays SSD force field
% 
%  7.3   HSITOGRAMS
% 
%      - E7_histogram1D                     creating a histogram
%      - E7_histogram1D_ext                 extended version of the above
% 
%  7.4  DISTANCES in MLPIR
% 
%      - E7_PETCT_MLPIR                     how to use distance.m in FAIR
%      - E7_PETCT_MLPIR_ext                 extended version of the above
% 
% ==============================================================================
% 
%  E8    REGULARIZATION
% 
%      - E8_Belastic 
%      - E8_checkOperations   test matrix based versus matrix free implementations
%      - E8_forces            compute response to forces
%      - E8_matrices          show and compare operators matrix free and matrix based
%      - E8_regularizationMB  how to use the matrix based regularizer 
%      - E8_regularizationMF  how to use the matrix free regularizer
%      - E8_setup             how to initialize parameters
% 
% ==============================================================================
% 
%  E9    NON-PARAMETRIC IMAGE REGISTRATION
% 
%  9.1   FIXED LEVEL
% 
%      - E9_HNSP_NPIR                   plain
%      - E9_HNSP_NPIR_pre               with pre-registration    
%      - E9_Hands_NPIR                  plain
%      - E9_Hands_NPIR_MI_mbElas_BFGS   using a lBFGS scheme          
%      - E9_Hands_affine                affine
%      - E9_Hands_NPIR_pre              with pre-registration  
%      - E9_Hands_NPIRmb_GN             using a Gauss-Newton scheme, matrix based
%      - E9_Hands_NPIRmf_GN             using a Gauss-Newton scheme, matrix free
%      - E9_Hands_NPIRmf_TR_nopre       using a Trust-Region scheme, matrix free
%      - E9_Hands_NPIRmf_TR_pcg         using a Trust-Region scheme, SOR preconditioned CG    
% 
%  9.2   MULTILEVEL
% 
%      - E9_HNSP_MLIR_SSD_mbCurv             SSD, matrix based curvature
%      - E9_HNSP_MLIR_SSD_mbElas             SSD, matrix based elastic
%      - E9_HNSP_MLIR_SSD_mfCurv             SSD, matrix free curvature
%      - E9_HNSP_MLIR_SSD_mfElas             SSD, matrix free elastic
%      - E9_HNSP_MLIR_TR                     SSD, using a Trust-Region method
%      - E9_Hands_MLIR_SSD_mbCurv            SSD, matrix based curvature
%      - E9_Hands_MLIR_SSD_mbElas            SSD, matrix based elastic
%      - E9_Hands_MLIR_SSD_mfCurv            SSD, matrix free curvature
%      - E9_Hands_MLIR_SSD_mfElas            SSD, matrix free elastic
%      - E9_MRIhead_MLIR_MI_mbElas           MI, matrix based elastic  
%      - E9_MRIhead_MLIR_NGF_mbElas          NGF, matrix based elastic
%      - E9_MRIhead_MLIR_SSD_mbElas          SSD, matrix based elastic
%      - E9_MRIhead_MLIRlBFGS_MI_mfElas      MI, matrix based elastic, using a lBFGS scheme    
%      - E9_MRIhead_MLIRlBFGS_NGF_elas       NGF, matrix based elastic, using a lBFGS scheme   
%      - E9_MRIhead_MLIRlBFGS_NGF_mbElas     NGF, matrix based elastic, using a lBFGS scheme
%      - E9_PETCT_MLIR_NGF_mbCurv            NGF, matrix based curvature 
%      - E9_PETCT_MLIR_NGF_mbElas            NGF, matrix based elastic 
%      - E9_PETCT_MLIR_NGF_mfElas            NGF, matrix free elastic  
%      - E9_PETCT_MLIRlBFGS_MI_mbCurv        MI, matrix based curvature, using a lBFGS scheme      
%      - E9_PETCT_MLIRlBFGS_MI_mbElas        MI, matrix based elastic, using a lBFGS scheme       
%      - E9_PETCT_MLIRlBFGS_NGF_mbCurv       NGF, matrix based curvature, using a lBFGS scheme       
%      - E9_PETCT_MLIRlBFGS_NGF_mbElas       NGF, matrix based elastic, using a lBFGS scheme       
% 
%  9.3   MULTISCALE
% 
%      - E9_Hands_MSIR             
% 
%  9.4   3D DATA
% 
%      - E9_3Dbrain_GN                      using a Gauss-Newton method 
%      - E9_3Dbrain_TR                      using a Trust-Region method
%      - E9_3Dbrain_lBFGS                   using a lBFGS        method
% 
% ==============================================================================
% 
%  U     AS SHOWN IN THE FAIR BOOK
% 
%      - book_affine2D            explaining affine transformation 
%      - book_kron3D              explaining Kronecker-products
%      - book_linearInter       explaining linear interpolation 
%      - book_rigid2D             explaining rigid transformation 
%      - book_splineInter1D       explaining spline interpolation 1D
%      - book_splineInter2D       explaining spline interpolation 2D
% 
% ==============================================================================
% 
%  S     SETUP DATA
% 
%      - get2Ddata                unified framework for setup*
%      - setup3DboxData           3D box  
%      - setup3DbrainData         3D brain    
%      - setupHNSPData            2D serial sectioning 
%      - setupHandData            2D x-ray 
%      - setupMRIData             2D MRI head
%      - setupMyData              generic example
%      - setupPETCTdata           2D PET / CT   
%      - setupUSData              2D ultra sound  
% 
% ==============================================================================
% 
%  X     MISCELLANEOUS           
% 
%      - FAIRdiary             unified framework for diaries/log-files
%      - pause             bypasses ''pause'' in the test-suite
% 
%      - tutWebsite            as shown on the website
% 
% ==============================================================================
% ============================================================================== 
