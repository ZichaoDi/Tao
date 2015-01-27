function g = getg_fdi(v)
%----------------------------------------------
% Compute the gradient g(v) 
% via finite-differencing with
% complex-valued perturbations.
%----------------------------------------------
% Usage: g = getg_fdi(v)
%----------------------------------------------
h = 1.0e-20;
g = v(:);
i = sqrt(-1);

for j = 1:length(v);
  vh    = v;
  vh(j) = vh(j) + i*h;
  fh    = getf(vh);
  g(j)  = imag(fh)/h;
end;