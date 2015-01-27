%----------------------------------------------------------------------
% Solve via traditional optimization
%----------------------------------------------------------------------
global current_n
global NF N        % NF counts # of function evals on each grid
global v_low v_up lambdatn
%----------------------------------------------------------------------
do_setup;
more off;
NF   = [0*N; 0*N; 0*N]; 
it   = 1;
fnl  = 0*v0;
n    = N(1);
global_setup(n);
current_n = n;
%----------------------------------------------------------------------
myfun = 'sfun';
e=cputime;
if (bounds);
   [v_low,v_up] = set_bounds(1,0);
  [v,F,G,ierror] = tnbc (v0,myfun,v_low,v_up);
   %[v,F,G,ierror] = tnbcm (v0,myfun,v_low,v_up,300);
else
   [v,F,G,ierror] = tn   (v0,myfun);
   %[v,F,G,ierror] = tnm   (v0,myfun,300);
end;
time=cputime-e;
disp(['time elpased is  ' , num2str(time)])
doplot(0,v);
report_results(N);
more on;