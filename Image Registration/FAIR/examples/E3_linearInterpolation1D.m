%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: linear interpolation in 1D, from the book
%
%==============================================================================

TD = [1,2,3,3]; m = length(T); omega = [0,4]; TD = reshape(TD,m,1);
xc = getCellCenteredGrid(omega,m);
Tc = book_linearInter1D(TD,omega,xc);
xf = linspace(omega(1)-1,omega(2)+1,201);
Tf = book_linearInter1D(TD,omega,xf);
clf; ph = plot(xc,TD,'b+',xc,Tc,'ro',xf,Tf,'k-');
