function [f,g]=sfun_XTM(W,M,MU_e,I0,Ltol,thetan,m,nTau,NumElement)
global SigMa_XTM LogScale eX eY
global Tik penalty XTMscale
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== M: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
W1=reshape(W,m(1),m(2),NumElement);
beta=1;
lambda=1e-3;
%%%%% =================== Attenuation Matrix at beam energy

MUe=reshape(MU_e(:,1,1),1,1,NumElement).*XTMscale;
MU=sum(W1.*repmat(MUe,[m(1),m(2),1]),3).*XTMscale;
%%%====================================
if(penalty)
    f=lambda*(norm(Tik*W))^2;
else
    f=0;
end
g=zeros(m(1),m(2),NumElement);
for n=1:length(thetan)
    sum_Tau=0;
    if(LogScale)
        Mt=-log(M(:,n)./I0);
        for i=1:nTau+1
            count=(nTau+1)*(n-1)+i;
            L=Ltol{n,i};
            if(~isempty(find(L,1)))
                Rdis=eX'*(MU.*L)*eY;
                sum_Tau=sum_Tau+beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g+2*beta*SigMa_XTM(count)*(Rdis-Mt(i)).*repmat(L,[1,1,NumElement]).*repmat(MUe,[m(1),m(2),1]);
            end
        end
    else
        Mt=M(:,n);
        for i=1:nTau+1
            count=(nTau+1)*(n-1)+i;
            L=Ltol{n,i};
            if(~isempty(find(L,1)))
                Rdis=I0*exp(-eX'*(MU.*L)*eY);%% Discrete case
                sum_Tau=sum_Tau+beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g-2*beta*SigMa_XTM(count)*Rdis*(Rdis-Mt(i)).*repmat(L,[1,1,NumElement]).*repmat(MUe,[m(1),m(2),1]);
            end
        end
    end
    f=f+sum_Tau;
end
g=g(:);
if(penalty)
    Wdiff=W(1:end-1)-W(2:end);
    g=g+lambda.*(-2.*[0;Wdiff]+2.*[Wdiff;0]);
    % g=g+lambda*2*Tik'*(Tik*W);
end





