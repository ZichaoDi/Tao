%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: multiscale spline interpolation in 1D
%
%==============================================================================

function varargout = E3_MS_splineInterpolation1D(regularizer)
if nargin == 0, regularizer = 'I';  end;                  % default value
dataT  = [0,2,2,2,1]'; m = length(dataT); omega = [0,m];  % initialize data
xc = linspace(omega(1)-1,omega(2)+1,101);                 % fine discretization 
Tc = @(T) book_splineInter1D(T,omega,xc);                 % spline interpolant
B  = spdiags(ones(m,1)*[1,4,1],[-1:1],m,m);               % spline basis
D  = spdiags(ones(m,1)*[-1,1],[0,1],m-1,m);               % derivative operation
M  = toeplitz([96,-54,0,6,zeros(1,m-4)]);                 % second derivative 
switch regularizer,                                       % initialize regularization   
  case 'I', W  = speye(m,m);                              % identity,   Tikhonov 
  case 'D', W  = D'*D;                                    % derivative, Tikhonov-Phillips
  case 'M', W  = M;                                       % bending operator
end;
c  = @(theta) (B'*B+theta*W)\(B'*dataT);  % coefficients as function in  theta
ph(1:2) = plot(getCellCenteredGrid(omega,m),dataT,'.k',xc,Tc(B\dataT),'k-');  hold on;
ph(3:5) = plot(xc,Tc(c(1)),'k-.',xc,Tc(c(19)),'k--',xc,Tc(c(100)),'k-');  hold off;
if nargout ~= 0; varargout = {ph}; end;