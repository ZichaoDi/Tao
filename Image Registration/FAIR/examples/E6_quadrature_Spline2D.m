%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: midpoint quadrature rule for 2D spline
%
%==============================================================================

clear, close all, help(mfilename)

omega = [0,6,0,8]; I = 36; T = zeros(6,8); T(3,4) = 1; h = []; Q = [];
psi = @(xc) splineInter(T,omega,xc);
for j=1:10,
  m    = 2^j*[1,1]; 
  h(j) = prod((omega(2:2:end)-omega(1:2:end))./m); 
  xc   = getCellCenteredGrid(omega,m); 
  Q(j) = h(j)*sum(psi(xc));
end;
figure(1); clf; p1=semilogx(h/h(1),Q+eps,'kx',h/h(1),Q,'k-');
figure(2); clf; p2=loglog(h/h(1),abs(I-Q)+eps,'kx',h/h(1),abs(I-Q),'k-');
