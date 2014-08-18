
global NF N 
global low up 
close all;
more off;
startup;
N=5;

NF = [0*N; 0*N; 0*N];
e=cputime;
low=0.05*ones(N^2,1);
up=1e6*ones(N^2,1);
% fctn=@(x)RosenBrook(x);
h=1/(N+1);
L = delsq(numgrid('S',N+2))/h^2;
xs=L\ones(N^2,1);
x0=xs+ones(N^2,1);
xinitial=x0;
fctn=@(x)Quadratic(x,L);
% foo(fctn,x0);
% return;
[xstar,f,g,ierror] = tnbc (x0,fctn,low,up);
doplot(xstar)
%%%====================================================== Report Result

t=cputime-e;
errTol=norm(xstar-xs)/norm(xinitial-xs);
fprintf('Time elapsed is %f, residule is %d\n',t,errTol);
%%%%%%%%%%%%%%%%%%%%=======================================================
