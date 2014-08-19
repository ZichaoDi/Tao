%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% For an extended documentation, see:
% Jan Modersitzki. FAIR: Flexible Algorithms for Image Registration, SIAM, 2009.
% http://www.siam.org/books/fa06/
% 
% Contents of FAIR data:
%
% This folder provides the data (see references at the end of the file) that come 
% with the toolbox and some convenient tools.
% 
% 2D Data:
%   hands-R.jpg, hands-T.jpg
%     xrays of hands                see reference [1] and setupHandData
%   HNSP-R.jpg,    HNSP-T.jpg     
%     images from serial sectioning see reference [2] and setupHNSPData
%   MRIhead-R.jpg, MRIhead-T.jpg  
%     MR slices from human brain    see reference [3] and setupMRIData
%   PET-CT-CT.jpg, PET-CT-PET.jpg
%     PET-CT study                  see reference [4] and setupPETCTdata
%   US Vibe Heldmann.jpg
%     ultrasound image of Vibe Heldmann
%   US.jpg
%     ultrasound image, http://www.mothercareultrasound.com/2d_ultrasound.jpg
%   EPIslice-R.jpg, EPIslice-T.jpg
%     2D slices of MR images, Institute of Clinical Radiology, Univesity
%     Hospital Muenster, Germany.
% 
% 3D Data:
%   brain3D.mat
%     MRI data of a human brain, data courtesy Ron Kikinis, Brigham & Women's Hospital, Boston, USA
%   createPhilipsKnees.mat
%     MRI data from a human knee, data courtesy Thomas Netsch, Philips Medical Solutions, Hamburg
%   mice3D.mat
%     PET data of a mouse heart, European Institute for Molecular Imaging (EIMI), University of Muenster, Germany
%   phantom3D.mat
%     MR data of a hardware phantom, Institute for Clinical Radiology, University Hospital Muenster, Germany
% 
% Mfiles in 'data':
% 	contents            - this file
%   checkSetupDataFile  - initializes data and parameters, if already available
% 	get2Ddata           - generic tool to generate the data structure needed
%   setup2Ddisc2CData   - an academical example: a disc and a C
%   setup2DEPIData      - 2D slices of human brain MR
%   setup2DGaussianData - an academical example: 2D Gaussian blobs
%   setup2DGaussianNoisyData - an academical example: noisy 2D Gaussian blobs
%   setup2DMPsquare     - an academical example: two squares with same mass
% 	setup3DboxData      - an academical example: a box in 3D
% 	setup3DbrainData    - human brain data
% 	setup3DkneeData     - initializing knee data
% 	setup3DmiceData     - initializing PET data of mouse heart
%   setup3DphantomData  - MR phantom data in 3D
% 	setupHandData       - xray's of human hands
% 	setupHNSPData       - images from histological serial sectioning
% 	setupMRIData        - MR slices from a human brain (T1/T2)
% 	setupMyData         - generic tool, explains how to get your data into FAIR
% 	setupPETCTData      - PET/CT images of human thorax's
% 	setupUSData         - ultrasound image (only one)
%
% data sources:
% [1]
%  @article{Amit1994,
%    author = {Yali Amit},
%     title = {A nonlinear variational problem for image matching},
%      year = {1994},
%   journal = {SIAM J. Sci. Comput.},
%    volume = {15},
%    number = {1},
%     pages = {207--224},
%  }
% [2] Oliver Schmitt, Institute of Anatomy, University of Rostock, Germany
% [3] http://www.bic.mni.mcgill.ca/brainweb/
% [4]
%  @article{ShekharEtAl2005,
%   author = {Raj Shekhar and Vivek Walimbe and Shanker Raja and Vladimir Zagrodsky
%             and Mangesh Kanvinde and Guiyun Wu and Bohdan Bybel},
%    title = {Automated 3-Dimensional Elastic Registration of Whole-Body {PET}
%             and {CT} from Separate or Combined Scanners},
%  journal = {J. of Nuclear Medicine},
%   volume = {46},
%   number = {9},
%     year = {2005},
%    pages = {1488--1496},
%   }
%
% % Copyright (c): Jan Modersitzki
% Version FAIR.2011
%==============================================================================
help(mfilename)