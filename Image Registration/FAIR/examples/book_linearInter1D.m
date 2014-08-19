%$  function Tc = linearInterpolation1D(T,omega,x);
%$ (c) Jan Modersitzki 2009/04/07, see FAIR.2 and FAIRcopyright.m.
%$ \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
function Tc = linearInterpolation1D(T,omega,x);

% get data size m, cell size h, dimension d, and number n of interpolation points
m  = length(T); 
h  = (omega(2:2:end)-omega(1:2:end))./m; 
d  = length(omega)/2; 
n  = length(x)/d;    
x  = reshape(x,n,d);
% map x from [h/2,omega-h/2] -> [1,m],
for i=1:d, x(:,i) = (x(:,i)-omega(2*i-1))/h(i) + 0.5; end;

Tc = zeros(n,1);                          % initialize output
Valid = @(j) (0<x(:,j) & x(:,j)<m(j)+1);  % determine indices of valid points
valid = find( Valid(1) );                
                                          % pad data to reduce cases
pad = 1; TP = zeros(m+2*pad,1); TP(pad+(1:m)) = reshape(T,m,1);
P = floor(x); x = x-P;                    % split x into integer/remainder
p = pad + P(valid); xi = x(valid);        % add the padding
Tc(valid) = TP(p).* (1-xi) + TP(p+1).*xi; % compute weighted sum
