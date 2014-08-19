%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: spline interpolation in 1D, from the book
%
%==============================================================================

dataT = [0,2,2,2,1]; m = length(dataT); omega = [0,m];  
xc = getCellCenteredGrid(omega,m);
B  = spdiags(ones(m,1)*[1,4,1],[-1:1],m,m);
T  = B\reshape(dataT,m,1);
xf = linspace(-1,6,101);             
Tf = book_splineInter1D(T,omega,xf);
figure(1); clf; ph = plot(xc,dataT,'.',xf,Tf,'g-','markersize',30); 
