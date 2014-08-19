%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR:  SSD versus discretization width (midpoint quadrature rule)
%
% - load data (setupHNSPData)
% - interpolate on various grids (splineInter)
% - compute SSD and plot it
%==============================================================================

clc, clear; help(mfilename);

setupHNSPData; clf; h = []; Q = [];
inter('reset','inter','splineInter');
[T,R] = inter('coefficients',dataT,dataR,omega,'out',0);
for j=1:10,
  m    = 2^j*[1,1]; 
  h(j) = prod((omega(2:2:end)-omega(1:2:end))./m); 
  xc   = getCellCenteredGrid(omega,m); 
  res  = inter(T,omega,xc) - inter(R,omega,xc);
  psi  = 0.5*h(j)*res'*res;
  Q(j) = psi;
end;
figure(1); clf; p1=semilogx(h/h(1),Q+eps,'kx',h/h(1),Q,'k-');
