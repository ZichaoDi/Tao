function [f,g]=sfun_XTM_com(M,MU,I0,Ltol,thetan,m,nTau,NumElement)

global Tik penalty
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== M: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
beta=1;
lambda=1e-3;
%%%====================================
MU=reshape(MU,m(1),m(2));
e=ones(m(1),1);
if(penalty)
f=lambda*(norm(Tik*W))^2;
else
    f=0;
end
g=zeros(m(1),m(2));
Jacob=g;
for n=1:length(thetan)
    Mt=-log(M(:,n)./I0);
%     Mt=M(:,n);
    sum_Tau=0;
    for i=1:nTau+1
        L=Ltol{n,i};
        if(~isempty(L))
            Rdis=e'*(MU.*L)*e;%
%             Rdis=I0*exp(-e'*(MU.*L)*e);%% Discrete case
            sum_Tau=sum_Tau+beta*(Rdis-Mt(i))^2;
            g=g+2*beta*(Rdis-Mt(i)).*L;%
            Jacob=Jacob+L;
            if(penalty)
%           g=g-2*beta*Rdis*(Rdis-Mt(i)).*repmat(full(L),[1,1,NumElement]).*repmat(MUe,[m(1),m(2),1]);

            end

        end
    end
    f=f+sum_Tau;
end

Jacob=reshape(Jacob,m(1),m(2));
condN=cond(Jacob);
[U,S,V]=svd(Jacob);
minS=min(diag(S));
maxS=max(diag(S));
fprintf('Condition %d, maximumS %d, minimumS %d \n',condN,maxS, minS);

g=g(:);
