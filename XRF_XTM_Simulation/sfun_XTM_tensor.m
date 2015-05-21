function [f,g,JJ]=sfun_XTM_tensor(DisR,MU,I0,Ltol,m,nTau)  %#codegen 
global Beta 
global numThetan 
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== DisR: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
MU=reshape(MU,1,1,m(1),m(2));
%%%====================================
Beta=1e0;%
I0=1.5e4;
r=I0*sum(sum(bsxfun(@times,MU,Ltol),3),4)-DisR';%+log(DisR'./I0);
g=sum(sum(2*Beta*bsxfun(@times,reshape(r,1,1,numThetan,nTau+1),permute(Ltol,[3,4,1,2])),3),4);
JJ=reshape(2*Beta*Ltol,numThetan*(nTau+1),prod(m));
f=Beta*sum(sum(r.^2,1),2);
% size(JJ)
% rank(JJ,1e-10)
% g=JJ'*r;
g=g(:);
% JJ=JJ(setdiff(1:numThetan*(nTau+1),EmptyBeam),:);
% f=f*prod(m);
% g=g(:)*prod(m);