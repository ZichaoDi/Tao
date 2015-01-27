function g = getg_fd(v)
%----------------------------------------------
% Compute the gradient g(v) 
% via finite-differencing with
% real-valued perturbations.
%----------------------------------------------
% Usage: g = getg_fd(v)
%----------------------------------------------
h = 1e-5;
g = v(:);
f = getf(v);

for j=1:length(v);
  vh    = v;
  vh(j) = vh(j) + h;
  fh    = getf(vh);
  g(j)  = (fh-f)/h;
end;