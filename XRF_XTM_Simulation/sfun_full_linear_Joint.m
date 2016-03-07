function [ft,gt,f_XRF,f_XTM]=sfun_full_linear(W,xrfData,DisR,MU_e,NumElement,m,nTau)
%%==== XRF objective function in least square form
%%==== Full dimension: nThetan * nTau * nv * ne * nE
%%==== Solve W when the W in exponential term is fixed 
global numChannel  numThetan 
global ConstSub 
mtol=prod(m);
% W=repmat(reshape(W,[1,1 mtol NumElement]),[numThetan,nTau+1 1 1 numChannel]);
W=reshape(W,[1,1 mtol NumElement]);
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,NumElement);
MUe_XTM=MUe*XTMscale;
MU_XTM=sum(bsxfun(@times,squeeze(W(1,:,:,1)),MUe),2)';
%%%%% ====================================================================
% XRF_v=squeeze(sum(sum(ConstSub.*W,3),4));
XRF_v=squeeze(sum(sum(bsxfun(@times,W,ConstSub),3),4));
%%======================================        
 ft=sum(sum(sum((XRF_v-xrfData).^2,3),2),1);
 g=2*sum(sum(sum(bsxfun(@times,ConstSub, ...
         reshape(XRF_v-xrfData,numThetan,nTau+1,1,1,numChannel)),1),2),5);
% [ft,g]=sfun_full_linearMexC(W,xrfData,NumElement,m,nTau,numChannel,numThetan,ConstSub);
gt=g(:);
