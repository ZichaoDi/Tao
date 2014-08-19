%==============================================================================
% (c) Jan Modersitzki 2011/04/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% For an extended documentation, see:
% Jan Modersitzki. FAIR: Flexible Algorithms for Image Registration, SIAM, 2009.
% http://www.siam.org/books/fa06/
% 
% Contents of FAIR TESTS
% 
%   contents         - this file
%   checkToolbox     - check a toolbox
%   getFairFiles     - returns all FAIR mfiles with patter [pattern] in folder [folder]
%   getFiles         - returns all mfiles with patter [pattern] in folder [folder]
%   parseDistance    - generic derivative test of the distance functions
%   Print            - check whether the function Print is still in the example
%   testDistance     - test distance measure tools
%   testDistancesMEX - test MEX distance measures
%   testDistancesOMP - test MEX distance measures and parallelization using OpenMP
%   testEnd          - 
%   testExamples     - test Examples and such schemes
%   testFAIR         - test FAIR package
%   testGeometryOMP  - test MEX geometry schemes and parallelization using OpenMP
%   testInter        - test interpolation schemes
%   testInterMEX     - test MEX interpolation schemes
%   testInterOMP     - test MEX interpolation schemes and parallelization using OpenMP
%   testKernel       - test kernel
%   testLM           - test landmark registration tools
%   testRegularizer  - test regualarization toolbox
%   testScriptFiles  - test script files
%   testSetupData    - test the setup of various examples
%   testStart        - 
%   testTrafo        - test transformation schemes
%==============================================================================
help(mfilename)