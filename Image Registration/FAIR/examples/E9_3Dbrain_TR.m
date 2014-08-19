% ===============================================================================
% Example for MLIR, MultiLevel Image Registration
% (c) Jan Modersitzki 2009/04/02, see FAIR.2 and FAIRcopyright.m.
% 
%   - data                 3d-brains, Omega=(0,20)x(0,10)x(0,20), level=3:6, m=[128,64,128]
%   - viewer               imgmontage
%   - interpolation        linearInter
%   - distance             SSD
%   - pre-registration     affine3D
%   - regularizer          mfElastic
%   - optimizer            Trust-Region
% ===============================================================================

setup3DbrainData;

% prepare the plot
FAIRplots('clear')
Dshow = @(T,R,omega,m) viewIP(abs(T-R),omega,m,'colormap',gray(256));
FAIRplots('set','Dshow',Dshow);

% finally: run the MultiLevel Non-Parametric Image Registration
[yc,wc,his] = MLIR(MLdata,'NPIR',@TrustRegion,'pfun',@MGsolver,'parametric',0);
