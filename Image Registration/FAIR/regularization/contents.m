% ==================================================================================
% (c) Jan Modersitzki 2010/12/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% FAIR Reguarization Toolbox
% 
% general purpose tools:
%   contents              - this file
%   regularizer           - the specific regularizer used in FAIR
%     initialize: regularizer('reset','regularizer','mbElastic',...
%                             'alpha',1e3,'mu',1,'lambda',0);
%     use:        [Sc,dS,d2S] = regularizer(yc-yRef,omega,m);
%
%  The regularizers for NPIR are:
%  ------------------------------
%  curvature            curvature regularizer (cell-centered grid)
%                           initialize: 'alpha'
%                           see, e.g., E10_2Ddisc2C_curvature
%  elastic              (linear) elastic regularizer (staggered grid)
%                           initialize: 'alpha', 'mu', 'lambda'
%                           see, e.g., E10_2Ddisc2C_elastic
%  hyperelastic         hyperelastic regularizer    (nodal grid)
%                           initialize: 'alphaLength', 'alphaArea', 'alphaVolume'
%                           see, e.g., E10_2Ddisc2C_hyperElastic
%
%  (linear) Differential operators are built in:
%
%  getCurvatureMatrix     - generates curvature regularizer matrix (cell-centered grid)
%  getElasticMatrixNodal  - generates elastic regularizer matrix (nodal grid)
%  getElasticMatrixStg    - generates elastic regularizer matrix (staggeres grid)
%  getGradientNodal       - generates gradient operator matrix (nodal grid)
%
%  For (nonlinear) hyperelastic regularization area and volume of tetrahedral partition 
%  are computed in
%  geometry               - contains areas and volumes functions
%                           used in hyperElastic
%  geometrymexC.h / .cpp  - MEX version of area and volume computation
%
% see also BigTutorialRegularizer
% ==================================================================================
help(mfilename)
