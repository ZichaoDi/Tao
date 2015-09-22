function [f,g,shift_yT]=sfun_Tw(W) 
global  SigMa_XTM numThetan NumElement LogScale 
global  L  m nTau DisR N MU_e XTMscale I0 
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== DisR: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
W1=reshape(W,m(1),m(2),NumElement);
%%%%% =================== Attenuation Matrix at beam energy
L=reshape(full(L),numThetan,nTau+1,m(1),m(2));
MUe=reshape(MU_e(:,1,1),1,1,NumElement).*XTMscale;
MU=sum(W1.*repmat(MUe,[m(1),m(2),1]),3).*XTMscale;
%%%====================================
eX=ones(m(1),1);
eY=ones(m(2),1);
%%%====================================
f=0;
g=zeros(m(1),m(2),NumElement);
shift_yT=zeros(numThetan,nTau+1);
for n=1:numThetan
    sum_Tau=0;
    if(LogScale)
        Mt=-log(DisR(:,n)/I0); 
        for i=1:nTau+1
            Lsub=reshape(L(n,i,:,:),m(1),m(2));
            count=(nTau+1)*(n-1)+i;
            if(~isempty(find(Lsub,1)))
                Rdis=eX'*(MU.*Lsub)*eY;
                sum_Tau=sum_Tau+1*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g+2*1*SigMa_XTM(count)*(Rdis-Mt(i)).*repmat(Lsub,[1,1,NumElement]).*repmat(MUe,[m(1),m(2),1]);
                shift_yT(n,i)=Rdis;
            end
        end
    else
        Mt=DisR(:,n);
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

