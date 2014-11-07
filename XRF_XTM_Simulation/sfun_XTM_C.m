function [f,g]=sfun_XTM_C(W,M,MU_e,I0,Ltol,thetan,m,nTau,NumElement)
global SigMa_XTM LogScale
global XTMscale
%%===== Reconstruction discrete objective based on matrix multiplication
%%===== L: intersection length matrix
%%===== M: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
W1=reshape(W,m(1),m(2),NumElement);
beta=1;
%%%%% =================== Attenuation Matrix at beam energy

MUe=reshape(MU_e(:,1,1),1,1,NumElement).*XTMscale;
MU=sum(W1.*repmat(MUe,[m(1),m(2),1]),3).*XTMscale;
%%%====================================
e1=ones(m(1),1);
g=zeros(prod(m)*NumElement,(nTau+1)*length(thetan));%zeros(m(1),m(2),NumElement);
sum_Tau=zeros(nTau+1,length(thetan));
for n=1:length(thetan)

    if(LogScale)
        Mt=-log(M(:,n)./I0);
        for i=1:nTau+1
            L=Ltol{n,i};
            if(~isempty(find(L,1)))
                Rdis=e1'*(MU.*L)*e1;
                sum_Tau(i,n)=Rdis-Mt(i);
%                 g=g+2*beta*SigMa_XTM(i)*(Rdis-Mt(i)).*repmat(L,[1,1,NumElement]).*repmat(MUe,[m(1),m(2),1]);
g(:,(nTau+1)*(n-1)+i)=reshape(repmat(L,[1,1,NumElement]).*repmat(MUe,[m(1),m(2),1]),prod(m)*NumElement,1);
            end
        end
    else
        Mt=M(:,n);
        for i=1:nTau+1
            L=Ltol{n,i};
            if(~isempty(find(L,1)))
                Rdis=I0*exp(-e1'*(MU.*L)*e1);%% Discrete case
                sum_Tau=sum_Tau+beta*SigMa_XTM(i)*(Rdis-Mt(i))^2;
                g=g-2*beta*SigMa_XTM(i)*Rdis*(Rdis-Mt(i)).*repmat(L,[1,1,NumElement]).*repmat(MUe,[m(1),m(2),1]);
            end
        end
    end
end
f=beta*sum(sum(sum_Tau'*SigMa_XTM*sum_Tau));

g=beta*g*reshape((SigMa_XTM*sum_Tau)',(nTau+1)*length(thetan),1);

