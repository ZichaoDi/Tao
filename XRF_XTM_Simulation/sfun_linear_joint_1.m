function [ft,gt,f_XRF,f_XTM]=sfun_linear_joint(W,xrfData,DisR,MU_e,L,NumElement,m,nTau)
%%==== XRF objective function in least square form
%%==== Full dimension: nThetan * nTau * nv * ne * nE
%%==== Solve W when the W in exponential term is fixed 
global numChannel  numThetan 
global ConstSub frame LogScale I0 s_a 
global SigMa_XTM SigMa_XRF
global Beta TempBeta penalty Tik lambda RealBeam
mtol=prod(m);
%%%%% =================== XRF 
XRF_v=ConstSub*W;
%%%%% =================== XRT and Attenuation Matrix at beam energy
MUe=squeeze(MU_e(:,1,1));
L=reshape(full(L),numThetan*(nTau+1),mtol);
MU_XTM=reshape(W,mtol,NumElement)*MUe;
if(s_a)
    Mt=reshape(DisR',numThetan*(nTau+1),1);
else
    Mt=-log(DisR'./I0)*1e0; 
    Mt=Mt(:)-min(Mt(:));
end
Rdis=L*MU_XTM;
xrfData=xrfData./1e1;
%%======================================        
 if(strcmp(frame,'EM'))
      thres=1;
      Rdis=Rdis+thres;
      Mt=Mt+thres;
      f_XTM=sum(-log(Rdis).*Mt+Rdis);
      g_XTM= L'*(-Mt./Rdis+1)*MUe';
      %%=======================================
      XRF_v=XRF_v+thres;
      xrfData=xrfData+thres;
      f_XRF=sum(-log(XRF_v).*xrfData+XRF_v);
      g=ConstSub'*(1-xrfData./XRF_v);
 elseif(strcmp(frame,'LS'))
      f_XTM=sum((Rdis-Mt).^2);
      g_XTM=2*L'*(Rdis-Mt)*MUe';
      %%=======================================
      f_XRF=sum((XRF_v-xrfData).^2);
      g=2*ConstSub'*(XRF_v-xrfData);
 elseif(strcmp(frame,'mix'))
      f_XTM=sum((Rdis-Mt).^2);
      g_XTM=2*L'*(Rdis-Mt)*MUe';
      %%=======================================
      XRF_v=XRF_v+1;
      xrfData=xrfData+1;
      f_XRF=sum(-log(XRF_v).*xrfData+XRF_v);
      g=ConstSub'*(1-xrfData./XRF_v);
 end
ft = TempBeta*f_XRF + Beta*f_XTM;
gt = TempBeta*g(:)  + Beta*g_XTM(:);
if(penalty)
    L1_norm=1;
    L2_norm=0;
    if(L1_norm)
        Reg=Tik*reshape(W,mtol,NumElement);
        ft=ft+lambda*(norm(Reg(:)))^2;
        gt=gt+lambda*2*reshape(Tik'*(Reg),mtol*NumElement,1);
    elseif(L2_norm)
        Reg=sum(abs(W(:)));
        ft=ft+lambda*Reg;
        gt=gt+lambda;
    end
end
