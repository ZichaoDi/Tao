function [f,g]=sfun_XTM(W,M,MU_e,I0,Ltol,thetan,m,nTau,NumElement)
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== M: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
W=reshape(W,m(1),m(2),NumElement);
%%%============== Rescale W to make unity contribution
DiscreteScale=0;
if(DiscreteScale)
for iN=1:NumElement
    W(:,:,iN)=W(:,:,iN)/MU_e(iN,1,1);
end
end
%%%====================================
beta=1;%1e4;
%%%%% =================== Attenuation Matrix at beam energy
MU=zeros(m);
MUe=reshape(MU_e(:,1,1),NumElement,1);
for i=1:m(1)
    for j=1:m(2)
        MU(i,j)=sum(reshape(W(i,j,:),NumElement,1).*MUe);
    end
end
MUe=reshape(MUe,1,1,NumElement);
e=ones(m(1),1);
f=0;
g=zeros(m(1),m(2),NumElement);

for n=1:length(thetan)
    Mt=-log(M(:,n)./I0);
    sum_Tau=0;
    for i=1:nTau+1
        L=Ltol{n,i};
        if(~isempty(L))
            Rdis=e'*(MU.*L)*e; %% Discrete case
            sum_Tau=sum_Tau+beta*(Rdis-Mt(i))^2;
            g=g+2*beta*(Rdis-Mt(i)).*repmat(full(L),[1,1,NumElement]).*repmat(MUe,[3,3,1]);
        end
    end
    f=f+sum_Tau;
end
g=g(:);