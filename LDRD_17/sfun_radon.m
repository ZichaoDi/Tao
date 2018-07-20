function [f,g,r]=sfun_radon(MU,Mt,Ltol) 
% [x,f,g,ierror] = tnbc (zeros(N^2,1),fctn,low(n_delta+1:end),up(n_delta+1:end)); % algo='TNbc';
global lambda frame N
%%===== Reconstruction discrete objective
%%===== Ltol: intersection length matrix
%%===== f: sum_i ||e^T(Ltol_i.*I)e-M_i||^2, i=1..theta
MU=MU(:);
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
penalty=0;
if(penalty)
    Tik=delsq(numgrid('S',N(1)+2)); 
    L1_norm=0;
    L2_norm=1;
    lambda=5e-7;
    if(L2_norm)
        Reg=Tik*MU;
        f=f+lambda*(sum(Reg.^2));
        g=g+lambda*2*Tik'*Tik*MU;
    elseif(L1_norm)
        Reg=sum(abs(W(:)));
        f=f+lambda*Reg;
        g=g+lambda;
    end
end


