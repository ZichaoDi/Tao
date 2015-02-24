function [f,g]=sfun_XTM_com(M,MU,I0,Ltol,thetan,m,nTau)
global SigMa_XTM LogScale Tol
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== M: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
beta=1;
MU=reshape(MU,m(1),m(2));
eX=ones(m(1),1);
eY=ones(m(2),1);
%%%====================================
f=0;
g=zeros(m(1),m(2));
for n=1:length(thetan)
    sum_Tau=0;
    if(LogScale)
        Mt=-log(M(:,n)./I0);
        for i=1:nTau+1
            count=(nTau+1)*(n-1)+i;
            L=reshape(Ltol(n,i,:),m(1),m(2))./Tol;
            if(~isempty(find(L,1)))
                Rdis=eX'*(MU.*L)*eY;
                sum_Tau=sum_Tau+beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;              
                g=g+2*beta*SigMa_XTM(count)*(Rdis-Mt(i)).*L;               
            end
        end
    else
        Mt=M(:,n);        
        for i=1:nTau+1
            count=(nTau+1)*(n-1)+i;
            L=reshape(Ltol(n,i,:),m(1),m(2))./Tol;
            if(~isempty(find(L,1)))
                Rdis=I0*exp(-eX'*(MU.*L)*eY);%% Discrete case
                sum_Tau=sum_Tau+beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g-2*beta*SigMa_XTM(count)*Rdis*(Rdis-Mt(i)).*L;
            end
        end
    end
    f=f+sum_Tau;
end
g=g(:);
