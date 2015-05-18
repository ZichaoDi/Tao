function [f,g,JJ]=sfun_XTM_com(DisR,MU,I0,Ltol,thetan,m,nTau)
global SigMa_XTM LogScale Beta EmptyBeam
global numThetan N
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== DisR: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
MU=reshape(MU,m(1),m(2));
eX=ones(m(1),1);
eY=ones(m(2),1);
%%%====================================
f=0;
% if(m(1)==N(1))
%     Beta=1;
% else
%     Beta=4;
% end
g=zeros(m(1),m(2));
JJ=zeros(numThetan*(nTau+1),prod(m));
r=zeros(numThetan*(nTau+1),1);
for n=1:numThetan
    sum_Tau=0;
    if(LogScale)
        Mt=-log(DisR(:,n)./I0);
        for i=1:nTau+1
            count=(nTau+1)*(n-1)+i;
            L=reshape(Ltol(n,i,:),m(1),m(2));
            if(~isempty(find(L,1)))
                Rdis=eX'*(MU.*L)*eY;
                sum_Tau=sum_Tau+Beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                r(i)=(Rdis-Mt(i));
                g=g+2*Beta*SigMa_XTM(count)*(Rdis-Mt(i)).*L; 
                JJ((n-1)*(nTau+1)+i,:)=reshape(2*Beta*SigMa_XTM(count)*L,1,prod(m));
            end
        end
    else
        Mt=DisR(:,n);        
        for i=1:nTau+1
            count=(nTau+1)*(n-1)+i;
            L=reshape(Ltol(n,i,:),m(1),m(2));
            if(~isempty(find(L,1)))
                Rdis=I0*exp(-eX'*(MU.*L)*eY);%% Discrete case
                sum_Tau=sum_Tau+Beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g-2*Beta*SigMa_XTM(count)*Rdis*(Rdis-Mt(i)).*L;
               JJ((n-1)*(nTau+1)+i,:)=reshape(2*Beta*SigMa_XTM(count)*Rdis*(Rdis-Mt(i)).*L,1,prod(m));

            end
            
        end
    end
    f=f+sum_Tau;
end
size(JJ)
rank(JJ,1e-10)
pause;
% g=JJ'*r;
g=g(:);
% JJ=JJ(setdiff(1:numThetan*(nTau+1),EmptyBeam),:);
% f=f*prod(m);
% g=g(:)*prod(m);