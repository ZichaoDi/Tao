%$  function Tc = splineInterpolation2D(T,omega,x);
%$ (c) Jan Modersitzki 2009/04/07, see FAIR.2 and FAIRcopyright.m.
%$ \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
function Tc = splineInterpolation2D(T,omega,x);

% get data size m, cell size h, dimension d, and number n of interpolation points
m  = size(T); 
h  = (omega(2:2:end)-omega(1:2:end))./m; 
d  = length(omega)/2; 
n  = length(x)/d;    
x  = reshape(x,n,d);
% map x from [h/2,omega-h/2] -> [1,m],
for i=1:d, x(:,i) = (x(:,i)-omega(2*i-1))/h(i) + 0.5; end;

Tc = zeros(n,1);                          % initialize output
Valid = @(j) (0<x(:,j) & x(:,j)<m(j)+1);  % determine indices of valid points
valid = find( Valid(1) & Valid(2) );                
if isempty(valid), return; end;
                                          % pad data to reduce cases
pad = 3; TP = zeros(m+2*pad); TP(pad+(1:m(1)),pad+(1:m(2))) = T;
P = floor(x); x = x-P;                    % split x into integer/remainder
p = @(j) P(valid,j); xi = @(j) x(valid,j);
i1 = 1; i2 = size(TP,1);                  % increments for linearized ordering
p  = (pad + p(1)) + i2*(pad + p(2) - 1);

for j1=-1:2,                              % Tc as weighted sum
  for j2=-1:2,
    Tc(valid) = Tc(valid) + TP(p+j1*i1+j2*i2).*b0(j1,xi(1)).*b0(j2,xi(2));
  end;
end;

function b = b0(j,xi)
switch j,
  case -1, b = (1-xi).^3;
  case 0,  b = -2*(1-xi).^3+xi.^3+6*(1-xi);
  case 1,  b = (1-xi).^3-2*xi.^3+6*xi;
  case 2,  b =  xi.^3;
end;

function b = db0(j,xi)
switch j,
  case -1, b = -3*(1-xi).^2;
  case 0,  b =  6*(1-xi).^2+3*xi.^2-6;
  case 1,  b = -3*(1-xi).^2-6*xi.^2+6;
  case 2,  b =  3*xi.^2;
end;
