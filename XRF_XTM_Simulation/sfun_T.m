function [f,g,JJ,shift_yT,eig_val]=sfun_T(MU)

global  SigMa_XTM numThetan
global LogScale Beta I0
global  L  m nTau DisR N
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== DisR: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
j=find(m(1)==N);
%  Beta=10^(8^(length(N)-j));%
MU=reshape(MU,m(1),m(2));
eX=ones(m(1),1);
eY=ones(m(2),1);
%%%====================================
f=0;
g=zeros(m(1),m(2));
JJ=zeros(numThetan*(nTau+1),prod(m));
shift_yT=zeros(numThetan,nTau+1);
for n=1:numThetan
    sum_Tau=0;
    if(LogScale)
        Mt=-log(DisR(:,n)./I0);
        for i=1:nTau+1
            count=(nTau+1)*(n-1)+i;
            Lsub=reshape(L(n,i,:),m(1),m(2));
            if(~isempty(find(Lsub,1)))
                Rdis=eX'*(MU.*Lsub)*eY;
                sum_Tau=sum_Tau+Beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g+2*Beta*SigMa_XTM(count)*(Rdis-Mt(i)).*Lsub;
               JJ((n-1)*(nTau+1)+i,:)=reshape(2*Beta*SigMa_XTM(count)*Lsub,1,prod(m));
                shift_yT(n,i)=Rdis;
            end
        end
    else
        Mt=DisR(:,n);
        for i=1:nTau+1
            count=(nTau+1)*(n-1)+i;
            Lsub=reshape(L(n,i,:),m(1),m(2));
            if(~isempty(find(L,1)))
                Rdis=I0*exp(-eX'*(MU.*Lsub)*eY);%% Discrete case
                sum_Tau=sum_Tau+Beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g-2*Beta*SigMa_XTM(count)*Rdis*(Rdis-Mt(i)).*Lsub;
                shift_yT(n,i)=Rdis;
            end
        end
    end
    f=f+sum_Tau;
end
 f=f*prod(m);
g=g(:)*prod(m);
