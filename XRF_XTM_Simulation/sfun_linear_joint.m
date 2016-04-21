function [ft,gt,f_XRF,f_XTM]=sfun_linear_joint(W,xrfData,DisR,MU_e,L,NumElement,m,nTau)
%%==== XRF objective function in least square form
%%==== Full dimension: nThetan * nTau * nv * ne * nE
%%==== Solve W when the W in exponential term is fixed 
global numChannel  numThetan 
global ConstSub EM LogScale I0 
global SigMa_XTM SigMa_XRF
global Beta TempBeta penalty Tik lambda
mtol=prod(m);
L=reshape(full(L),numThetan,nTau+1,mtol);
%%%%% =================== Attenuation Matrix at beam energy
W=reshape(W,[1,1 mtol NumElement]);
% MUe=ones(size(squeeze(MU_e(:,1,1))));
MUe=squeeze(MU_e(:,1,1));
MU_XTM=squeeze(W)*MUe;
SigMa_XTM=1./DisR';
SigMa_XRF=1./xrfData;
SigMa_XTM=ones(size(SigMa_XTM));
SigMa_XRF=ones(size(SigMa_XRF));
Mt=-log(DisR'./I0);%DisR'*1e0;% 
Rdis=sum(bsxfun(@times,reshape(MU_XTM,[1,1,mtol]),L),3); %% Discrete case
%%%%% ====================================================================
XRF_v=squeeze(sum(sum(bsxfun(@times,W,ConstSub),3),4));
%%======================================        
if(EM)
     Rdis = Rdis +1;
     Mt = Mt + 1;
     f_XTM=sum(sum((-log(Rdis).*Mt+Rdis),2),1);
     g_XTM =squeeze(sum(sum(bsxfun(@times,reshape(bsxfun(@times,1-Mt./Rdis,L),1,numThetan,nTau+1,mtol),MUe),2),3))';
     XRF_v=XRF_v+1;
     xrfData=xrfData+1;
     f_XRF=sum(sum(sum((-log(XRF_v).*xrfData+XRF_v),3),2),1);
     g=1*sum(sum(sum(bsxfun(@times,ConstSub, ...
         reshape(1-xrfData./XRF_v,numThetan,nTau+1,1,1,numChannel)),1),2),5);
else
     f_XRF=sum(sum(sum((XRF_v-xrfData).^2.*SigMa_XRF,3),2),1);
     g=2*sum(sum(sum(bsxfun(@times,ConstSub, ...
         reshape((XRF_v-xrfData).*SigMa_XRF,numThetan,nTau+1,1,1,numChannel)),1),2),5);
     f_XTM=sum((Rdis(:)-Mt(:)).^2.*SigMa_XTM(:));
     g_XTM =2*squeeze(sum(sum(bsxfun(@times,reshape(bsxfun(@times,(Rdis-Mt).*SigMa_XTM,L),1,numThetan,nTau+1,mtol),MUe),2),3))';
end
ft = TempBeta*f_XRF+Beta*f_XTM;
gt=Beta*g_XTM(:)+TempBeta*g(:);
if(penalty)
    L1_norm=1;
    L2_norm=0;
    if(L1_norm)
        Reg=Tik*squeeze(W);
        ft=ft+lambda*(norm(Reg(:)))^2;
        gt=gt+lambda*2*reshape(Tik'*(Reg),mtol*NumElement,1);
    elseif(L2_norm)
        Reg=sum(abs(W(:)));
        ft=ft+lambda*Reg;
        gt=gt+lambda;
    end
end
