function [f,g,r]=sfun_radon(MU,Mt,I0,Ltol) 
global frame
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
if(strcmp(frame,'EM'))
    thres=1;
    Rdis=Ltol*MU+thres;
    Mt=Mt(:)+thres;
    f=sum(-log(Rdis).*Mt+Rdis);
    g= Ltol'*(-Mt./Rdis+1);
elseif(strcmp(frame,'LS'))
    r=Ltol*MU-Mt(:);
    g=Ltol'*r;
    f=1/2*sum(r.^2,1);
end

