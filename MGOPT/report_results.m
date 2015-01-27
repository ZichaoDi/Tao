function report_results(N)
%-------------------------------------------
% Print the grid sizes, as well as the
% # of function evaluations on each grid.
%
% The grid sizes are in N, the totals in NF.
%-------------------------------------------
% Usage: report_results(N)
%-------------------------------------------
global NF
%-------------------------------------------
fprintf('\n')
drawline(N);
fprintf('Optimization costs per grid\n')

drawline(N);

fprintf('N:     ');
fprintf(' %5i',N);
fprintf('\n');

drawline(N);

for j=1:3;
  if (j==1); fprintf('NF(it):'); end;
  if (j==2); fprintf('NF(nf):'); end;
  if (j==3); fprintf('NF(cg):'); end;
  fprintf(' %5i',NF(j,:));
  fprintf('\n');
end;

drawline(N);
fprintf('\n');
%===========================================
function drawline(N)

fprintf('-------');
for i=1:length(N);
   fprintf('------');
end;
fprintf('\n');