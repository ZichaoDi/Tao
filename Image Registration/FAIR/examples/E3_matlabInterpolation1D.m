%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: using MATLAB's linear interpolation in 1D 
%
%==============================================================================

omega = [0,4]; T = [1,2,3,3]'; m = length(T);
xc = getCellCenteredGrid(omega,m);
Tc = linearInterMatlab(T,omega,xc);
xf = linspace(omega(1)-1,omega(2)+1,201);
Tf = linearInterMatlab(T,omega,xf);
clf; ph = plot(xc,0*Tc,'b+',xc,Tc,'ro',xf,Tf,'k-');

fctn = @(x) linearInterMatlab(T,omega,x);
fig = checkDerivative(fctn,xc+1e-2); pause(2); close(fig);
