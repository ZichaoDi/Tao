function [f,g, JJ]=sfun_XTM_tensor(Mt,MU,I0,Ltol,m,nTau)  %#codegen 
global numThetan synthetic  Tol
global pert frame
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== DisR: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
frame='EM';
if(strcmp(frame,'EM'))
    thres=1;
    Rdis=Ltol*MU+thres;
    Mt=Mt(:)+thres;
    f=sum(-log(Rdis).*Mt+Rdis);
    g= Ltol'*(-Mt./Rdis+1);
else
    r=Ltol*MU-pert-Mt(:);
    g=Ltol'*r;
    f=1/2*sum(r.^2,1);
end
JJ = reshape(Ltol,numThetan*(nTau+1),prod(m));

