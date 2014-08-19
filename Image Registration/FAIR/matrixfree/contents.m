%==============================================================================
% (c) Jan Modersitzki 2011/04/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% For an extended documentation, see:
% Jan Modersitzki. FAIR: Flexible Algorithms for Image Registration, SIAM, 2009.
% http://www.siam.org/books/fa06/
% 
% Contents of FAIR MATRIXFREE
%
%   MGsolver - calls MG solver
%   mfvcyce  - MG solver
%   mfPu.m   - matrix free prolongation operator
%   mfJacobi - MG smoother
%   mfAy     - MG matrix vector operation (M+alpha*B'*B)
%   mfBy     - matrix free B*y for discrete regularizer
%==============================================================================
help(mfilename)
