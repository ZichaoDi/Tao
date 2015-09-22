%----------------------------------------------------------------------
% Solve with successive refinement:
% - optimize on coarsest grid, update, and repeat up to finest grid
%----------------------------------------------------------------------
global current_n v_low v_up
global NF N maxiter % NF counts # of function evals on each grid
%----------------------------------------------------------------------
% do_setup;
more off;
NF   = [0*N; 0*N; 0*N];
it   = 1;
fnl  = 0*v0;
current_n = N(end);
global_setup;
current_n = N(1);
%----------------------------------------------------------------------
v0=xinitial;
for j=2:length(N)
  v0        = downdate(v0,0);
  vH0=v0;
  current_n = N(j);
end;
for j=length(N):-1:2;
   fprintf('=======================\n');
   fprintf('Optimizing for n = %4i\n',N(j));
   fprintf('=======================\n');
   nit_solve=100;
   if (bounds);
      [v_low,v_up] = set_bounds(j);
      nit = nit_solve; [v0,F,G,ierror] = tnbcm (v0,'sfun',v_low,v_up,nit);
   else
      nit = nit_solve; [v0,F,G,ierror] = tnm   (v0,'sfun',nit);
   end;
   current_n = N(j-1);
   vHs=v0;
   v0 = update(v0,0);
   vupdate=v0;
   global_setup;
end;
fprintf('=======================\n');
fprintf('Optimizing for n = %4i\n',N(1));
fprintf('=======================\n');
maxiter=200;
if (bounds);
  [v_low,v_up] = set_bounds(1);
  [v,F,G,ierror] = tnbc (v0,'sfun',v_low,v_up);
else
  [v,F,G,ierror] = tn   (v0,'sfun');
end;

doplot(it,v,W_level);
report_results(N);
more on;
