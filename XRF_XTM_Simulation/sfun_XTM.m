function [f,g]=sfun_XTM(W,DisR,MU_e,I0,L,m,NumElement)
global Tik penalty frame 
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== M: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
mtol=prod(m);
lambda=1e-3;
%%%%% =================== Attenuation Matrix at beam energy
MUe=squeeze(MU_e(:,1,1));
MU_XTM=reshape(W,mtol,NumElement)*MUe;
Rdis=L*MU_XTM;
%%%====================================
if(penalty)
    f=lambda*(norm(Tik*W))^2;
else
    f=0;
end
if(strcmp(frame,'EM'))
    Mt=Mt+1;
    Rdis=Rdis+1;
     f_XTM=sum(-log(Rdis).*Mt+Rdis);
     g_XTM= L'*(-Mt./Rdis+1)*MUe';
elseif(strcmp(frame,'LS'))
     f_XTM=sum((Rdis-Mt).^2);
     g_XTM=2*L'*(Rdis-Mt)*MUe';
end
g=g_XTM(:);
f=f+f_XTM;
if(penalty)
    Wdiff=W(1:end-1)-W(2:end);
    g=g+lambda.*(-2.*[0;Wdiff]+2.*[Wdiff;0]);
    % g=g+lambda*2*Tik'*(Tik*W);
end





