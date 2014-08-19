%==============================================================================
% (c) Jan Modersitzki 2011/04/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% For an extended documentation, see:
% Jan Modersitzki. FAIR: Flexible Algorithms for Image Registration, SIAM, 2009.
% http://www.siam.org/books/fa06/
% 
% Contents of FAIR VIEWER TOOLBOX
% 
% General Viewer Tools:
%   contents               - this file
%   FAIRfigure             - opens a non gray FAIR figure at standard position
%   FAIRplots              - generates FAIR plots
%   FAIRposition           - computes the default position of a figure in FAIR
%   FAIRprint              - place to call a sophisticated printing tool
%   plotGrid               - Plot a d-dimensional grid
%   plotIterationHistory   - plots iteration history
%   plotMLIterationHistory - plots multi-level iteration history
%   rgb2gray               - converts RGB images to grayscale
%   showResults            - visualizes registration results
%   viewImage              - Main function for image visualization
% 
% 2D Viewer:
%   imgOverlay             - overlay of two images
%   overlayImage2D         - 2D image viewer for overlaying images R and T
%   viewImage2D            - THE 2D image viewer
%   viewImage2Dsc          - 2D image viewer
% 
% 3D Viewer:
%   imgmontage             - visualizes a 3D image as montage
%   scrollGrids            - Scrollable viewer for one or multiple 3D grid(s)
%   scrollImages           - Scrollable viewer for one or multiple 3D volume(s)
%   scrollOverlay          - Scrollable overlay function for two 3D volumes
%   scrollVolume           - scroll and compare 3D images
%   viewIP                 - visualizes the intensity projections of a 3D image
%   viewSlices             - visualizes a 3D image as slides
%   volView                - 3D image viewer
%==============================================================================
help(mfilename)