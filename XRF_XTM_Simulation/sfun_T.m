function [f,g,JJ,shift_yT]=sfun_T(MU)

global  SigMa_XTM numThetan
global LogScale  I0 
global  L  m nTau DisR N
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== DisR: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
scale=prod(m);
DR=DisR';
%%%====================================
precond=0;
if(precond)
    r=L'*(L*MU+log(DR(:)/I0));
    g=L'*L*r;
    JJ = L'*L*scale;
else 
    r=L*MU+log(DR(:)/I0);
    g=r'*L;
    JJ = reshape(L,numThetan*(nTau+1),prod(m))*scale;
end
f=1/2*sum(sum(r.^2,1),2)*scale^2;
g=g(:)*scale^2;
shift_yT=reshape(L*MU,numThetan,nTau+1);
%%%====================================
% j=find(m(1)==N);





