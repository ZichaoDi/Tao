function [f,g]=sfun_XTM_com(M,MU,I0,Ltol,thetan,m,nTau)
global SigMa_XTM LogScale
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== M: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
beta=1;
MU=reshape(MU,m(1),m(2));
e=ones(m(1),1);
f=0;
g=zeros(m(1),m(2));
for n=1:length(thetan)
    sum_Tau=0;
    if(LogScale)
        Mt=-log(M(:,n)./I0);
        for i=1:nTau+1
            L=Ltol{n,i};
            if(~isempty(find(L,1)))
                Rdis=e'*(MU.*L)*e;
                sum_Tau=sum_Tau+beta*SigMa_XTM(i)*(Rdis-Mt(i))^2;
                g=g+2*beta*SigMa_XTM(i)*(Rdis-Mt(i)).*L;               
            end
        end
    else
        Mt=M(:,n);        
        for i=1:nTau+1
            L=Ltol{n,i};
            if(~isempty(find(L,1)))
                Rdis=I0*exp(-e'*(MU.*L)*e);%% Discrete case
                sum_Tau=sum_Tau+beta*SigMa_XTM(i)*(Rdis-Mt(i))^2;
                g=g-2*beta*SigMa_XTM(i)*Rdis*(Rdis-Mt(i)).*repmat(full(L),[1,1,NumElement]).*repmat(MUe,[m(1),m(2),1]);
            end
        end
    end
    f=f+sum_Tau;
end
g=g(:);
