function [f,g, JJ]=sfun_XTM_tensor(DisR,MU,I0,Ltol,m,nTau)  %#codegen 
global numThetan synthetic  Tol
global pert
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== DisR: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
scale=1e6;
DR=-log(DisR'./I0);
r=Ltol*MU-pert-DR(:);
g=Ltol'*r;
JJ = reshape(Ltol,numThetan*(nTau+1),prod(m));
f=1/2*sum(r.^2,1)*scale;
g=g*scale;
