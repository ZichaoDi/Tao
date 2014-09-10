% [W0,XRF,RMlocal,I,M,Energy,EnergyChannel,Ltol,GlobalInd,thetan,m,nTau]=XRFreconCompond;
% save('data/Sphere4_90XRF.mat','W0','XRF','RMlocal','I','M','Energy','EnergyChannel','Ltol','GlobalInd','thetan','m','nTau');
% return;
more off;
clc;
load Sphere4_90XRF;
fctn=@(W)sfun_XRF(W,XRF,RMlocal,I,M,Energy,EnergyChannel,Ltol,GlobalInd,thetan,m,nTau);
x0=0.25*rand(prod(m)*size(M,1),1);%I0(:)+0*rand(size(I0(:)));%
%  x0=0.005*rand(prod(m)*size(M,1),1).*W0(:);
%  foo(fctn,x0);
% retun;
options = optimset('GradObj','on','Display','iter');
[x,fval] = fmincon(fctn,x0,[],[],[],[],zeros(size(x0)),ones(size(x0)),[],options);
err=norm(x-W0(:));
fprintf('residule is %d \n',err);

% OPT on Sphere4_90XRF:
%                                 Norm of      First-order 
%  Iteration        f(x)          step          optimality   CG-iterations
%      0        8.13795e+06                      9.47e+06                
%      1             2009.8        2.74326       4.76e+04           2
%      2           0.646193        1.23884       2.02e+03           5
%      3          0.0050833     0.00558058           97.7           3
%      4        5.38529e-05     0.00268827           13.2           3
%      5        4.90299e-08    6.95347e-05          0.256           5
%      6        3.66807e-11    1.17277e-05         0.0131           5
% 
% Local minimum possible.
% 
% fminunc stopped because the final change in function value relative to 
% its initial value is less than the default value of the function tolerance.
% 
% <stopping criteria details>
% 
% residule is 4.640390e-08 
