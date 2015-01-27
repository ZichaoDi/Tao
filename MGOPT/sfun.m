function [f, g] = sfun (v)
%--------------------------------------------------------------
% compute objective function and gradient
%--------------------------------------------------------------
global grad_type
%--------------------------------------------------------------
f = getf(v);
if (strcmp(grad_type,'adj'));
  g = getg_adj(v);        % via adjoint, or exact
end;
if (strcmp(grad_type,'fdr'));
  g = getg_fd (v);        % via finite-difference [real]
end;
if (strcmp(grad_type,'fdi'));
  g = getg_fdi(v);        % via finite-difference [complex]
end