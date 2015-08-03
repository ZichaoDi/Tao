function [f,g]=sfun_XTM(W,M,MU_e,I0,L,thetan,m,nTau,NumElement)
global SigMa_XTM LogScale numThetan 
global Tik penalty XTMscale 
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== M: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
W1=reshape(W,m(1),m(2),NumElement);
lambda=1e-3;
%%%%% =================== Attenuation Matrix at beam energy

MUe=reshape(MU_e(:,1,1),1,1,NumElement).*XTMscale;
MU=sum(W1.*repmat(MUe,[m(1),m(2),1]),3).*XTMscale;
%%%====================================
eX=ones(m(1),1);
eY=ones(m(2),1);
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
        Mt=-log(M(:,n)).*1e-7;%M(:,n)./1e7;% 
        for i=1:nTau+1
            Lsub=full(reshape(L(sub2ind([numThetan,nTau+1],n,i),:),m(1),m(2)));
            count=(nTau+1)*(n-1)+i;
            if(~isempty(find(Lsub,1)))
                Rdis=eX'*(MU.*Lsub)*eY;
                sum_Tau=sum_Tau+1*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g+2*1*SigMa_XTM(count)*(Rdis-Mt(i)).*repmat(Lsub,[1,1,NumElement]).*repmat(MUe,[m(1),m(2),1]);
            end
        end
    else
        Mt=M(:,n);
        for i=1:nTau+1
            Lsub=reshape(L(n,i,:),m(1),m(2));
            count=(nTau+1)*(n-1)+i;
            if(~isempty(find(Lsub,1)))
                Rdis=I0*exp(-eX'*(MU.*Lsub)*eY);%% Discrete case
                sum_Tau=sum_Tau+1*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g-2*1*SigMa_XTM(count)*Rdis*(Rdis-Mt(i)).*repmat(Lsub,[1,1,NumElement]).*repmat(MUe,[m(1),m(2),1]);
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





