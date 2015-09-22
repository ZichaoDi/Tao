function [f,g, JJ]=sfun_XTM_tensor(DisR,MU,I0,Ltol,m,nTau)  %#codegen 
global numThetan 
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== DisR: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
scale=prod(m);
DR=DisR';
%%%====================================
precond=0;
if(precond)
    r=Ltol'*(Ltol*MU+log(DR(:)/I0));
    g=Ltol'*Ltol*r;
    JJ = Ltol'*Ltol*scale;
else 
    r=Ltol*MU+log(DR(:)/I0);
    g=r'*Ltol;
    JJ = reshape(Ltol,numThetan*(nTau+1),prod(m))*scale;
end
f=1/2*sum(sum(r.^2,1),2)*scale^2;
g=g(:)*scale^2;
