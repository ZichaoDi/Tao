% (c) Jan Modersitzki 2010-08-11, see FAIR.2 and FAIRcopyright.m.
% \url{http://www.mic.uni-luebeck.de/people/jan-modersitzki/software/fair.html}
%
% FAIR KERNEL
%
%  contents              this file
%  checkDerivative       enables a test of derivatives
%  options               tool for handling of persistent variables
%  reportStatus          reports status of current FAIR setting 
% 
%  optim                 driver for optimization
%  setOptPara            setup of optimization parameter
%  Armijo                Armijo's line search, THE line search
%  ArmijoBacktrack       Armijo's line search with backtracking
%  ArmijoBacktrackPIR    Armijo's line search with backtracking for PIR
%  GaussNewton           Generalized Gauss-Newton scheme, THE optimizer
%  SteepestDescent       Steepest descent, downhill, ... 
%  TrustRegion           Trust region optimizer
%  lBFGS                 limited BFGS optimizer
%  xSquare               test function for optimization
%
%  MLIR                  Multi-level Image Registration THE scheme
%  MLPIR                 Multi-level Parametric Image Registration 
%  NPIRBFGSobjFctn       objectiv fcuntion for non-parametric image registration and BFGS
%  NPIRobjFctn           objectiv fcuntion for non-parametric image registration
%  PIRBFGSobjFctn        objectiv fcuntion for parametric image registration and BFGS
%  PIRobjFctn            objectiv fcuntion for parametric image registration
%  center                projection operator from any grid to cell centered
%  grid2grid             projection operator, interpolates from any grid to any
%  nodal2center          maps nodal to cell centered grid
%  stg2center            maps staggered to cell centered grid
%
%  getCellCenteredGrid   computes a cell-centered grid of size m of domain omega
%                        replaced previous version getCenteredGrid
%  getCenteredGrid       old version of getCellCenteredGrid
%  getNodalGrid          computes a nodal grid of size m of domain omega
%  getStaggeredGrid      computes a staggered grid of size m of domain omega
%  getMultilevel         provides a multi-level representation of the data
help(mfilename)