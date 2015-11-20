function [f,g,JJ,shift_yT]=sfun_T(MU)

global  SigMa_XTM numThetan
global LogScale  I0 
global  L  m nTau DisR N Lap
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== DisR: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
j=find(m(1)==N);
scale=prod(m-1);
DR=DisR';
%%%====================================
rng('default');
noise=0*1e8*rand(size(DR(:)));
r=L*MU+log((DR(:)+noise)/I0);
g=r'*L;
JJ = reshape(L,numThetan*(nTau+1),prod(m))*scale;
Tik=Lap{j};
beta=0*1e-7;
alpha=1;
f=1/2*sum(sum(r.^2,1),2)*scale^2;
g=g(:)*scale^2;
%%===1-norm regularizer
%f=alpha*f+beta*norm(Tik*MU,1);
%g=alpha*g+2*beta*Tik'*sign(Tik*MU);
%%===2-norm regularizer
f=alpha*f+beta*MU'*Tik'*Tik*MU;
g=alpha*g+2*beta*Tik'*Tik*MU;
shift_yT=reshape(L*MU,numThetan,nTau+1);
%%%====================================








