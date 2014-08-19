%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: midpoint quadrature rule for 1D spline
%
%==============================================================================

clear, close all, help(mfilename)

omega = [0,6];   I = 1;  psi = @(x) spline1D(3,x); h = []; Q = [];
for j=1:10,
  m    = 2^j; 
  h(j) = diff(omega)/m; 
  xc   = getCellCenteredGrid(omega,m); 
  Q(j) = h(j)*sum(psi(xc));
end;
figure(1); clf; p1=semilogx(h/h(1),Q+eps,'kx',h/h(1),Q,'k-');
figure(2); clf; p2=loglog(h/h(1),abs(I-Q)+eps,'kx',h/h(1),abs(I-Q),'k-');
