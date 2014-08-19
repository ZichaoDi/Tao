%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% For an extended documentation, see:
% Jan Modersitzki. FAIR: Flexible Algorithms for Image Registration, SIAM, 2009.
% http://www.siam.org/books/fa06/
% 
% Contents of FAIR Transformation Toolbox
%
%   contents    - this file
%   trafo       - the specific transformationin FAIR
%     initialize: trafo('reset','trafo','rigid2D');
%     usage:      w0 = trafo('w0'); % returns parameterization of identity
%                 [yc,dy] = trafo(wc,xc;
%
% affine2D                     - affine linear transformation for 2D
% affine2Dsparse               - affine linear transformation for 2D (memory efficient)
% affine3D                     - affine linear transformation for 3D
% affine3Dsparse               - affine linear transformation for 3D (memory efficient)
% matVecQw                     - matrix-vector multiplication for efficient versions
% rigid2D                      - rigid  transformation for 2D
% rigid2Dsparse                - rigid  transformation for 2D (memory efficient)
% rigid3D                      - rigid  transformation for 3D
% rigid3Dsparse                - rigid  transformation for 3D (memory efficient)
% rotation2D                   - rotation for 2D
% splineTransformation2D       - spline transformation for 2D
% splineTransformation2Dsparse - spline transformation for 2D (memory efficient)
% splineTransformation3Dsparse - spline transformation for 3D (memory efficient)
% translation2D                - translation for 2D
% translation3D                - translation for 3D
% tensorProdC.c                - C implementation of kron(Q3,Q2,Q1) * w
%                                see splineTransformation3Dsparse.m
%  
%  see also E6_HNSP_PIR_SSD_rotation2D_level4 and BigTutorialTrafo
%==============================================================================
help(mfilename)