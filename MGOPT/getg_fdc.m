function g = getg_fdc(v)
%----------------------------------------------
% Compute the gradient g(v) 
% via central-differencing with
% real-valued perturbations.
%----------------------------------------------
% Usage: g = getg_fd(v)
%----------------------------------------------
h = 1e-6;
g = v(:);
f = getf(v);

for j=1:length(v);
  vh    = v;
  vh(j) = vh(j) + h;
  fh    = getf(vh);
  vh    = v;
  vh(j) = vh(j) - h;
  fh1   = getf(vh);
  g(j)  = (fh-fh1)/(2*h);
end;
