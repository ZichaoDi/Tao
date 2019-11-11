function [f,g,r]=sfun_radon(MU,Mt,Ltol) 
% [x,f,g,ierror] = tnbc (zeros(N^2,1),fctn,low(n_delta+1:end),up(n_delta+1:end)); % algo='TNbc';
global f1 f2 reg lambda frame N
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
penalty=1;
if(penalty)
    lambda=1e-4;
    if(strcmp(reg,'L2'))
        Tik=delsq(numgrid('S',N(1)+2)); 
        Reg=Tik*MU;
        f2=sum(Reg.^2);
        f1=f;
        f=f+lambda*(sum(Reg.^2));
        g=g+lambda*2*Tik'*Tik*MU;
    elseif(strcmp(reg,'L1'))
        Reg=sum(abs(MU(:)));
        f=f+lambda*Reg;
        g=g+lambda;
    elseif(strcmp(reg,'TV'))
        epsilon=1e-3;
        f1=f;
        % Tik=diag([-1;ones(N-2,1);1])+diag(1/2*ones(N-1,1),1)+diag(-1/2*ones(N-1,1),-1);
        [Dx,Dy]=VecDiv(N);
        ux=Dx*MU;
        uy=Dy*MU;
        sum_Ux=sum(ux.^2);
        sum_Uy=sum(uy.^2);

        f2 = sqrt(sum_Ux + sum_Uy + epsilon^2);
        f = f+lambda*f2;
        shiftg = (Dx'*ux+Dy'*uy)./f2;
        g = g+lambda*shiftg(:);
    end
end


