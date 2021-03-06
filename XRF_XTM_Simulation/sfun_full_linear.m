function [ft,gt]=sfun_full_linear(W,xrfData,NumElement,m,nTau)
%%==== XRF objective function in least square form
%%==== Full dimension: nThetan * nTau * nv * ne * nE
%%==== Solve W when the W in exponential term is fixed 
global numChannel  numThetan 
global ConstSub EM 
mtol=prod(m);
% W=repmat(reshape(W,[1,1 mtol NumElement]),[numThetan,nTau+1 1 1 numChannel]);
W=reshape(W,[1,1 mtol NumElement]);
%%%%% ====================================================================
% XRF_v=squeeze(sum(sum(ConstSub.*W,3),4));
XRF_v=squeeze(sum(sum(bsxfun(@times,W,ConstSub),3),4));
%%======================================        
if(EM)
     XRF_v=XRF_v+1;
     xrfData=xrfData+1;
     ft=sum(sum(sum((-log(XRF_v).*xrfData+XRF_v),3),2),1);
     g=1*sum(sum(sum(bsxfun(@times,ConstSub, ...
         reshape(1-xrfData./XRF_v,numThetan,nTau+1,1,1,numChannel)),1),2),5);
else
     ft=sum(sum(sum((XRF_v-xrfData).^2,3),2),1);
     g=2*sum(sum(sum(bsxfun(@times,ConstSub, ...
         reshape(XRF_v-xrfData,numThetan,nTau+1,1,1,numChannel)),1),2),5);
end
% [ft,g]=sfun_full_linearMexC(W,xrfData,NumElement,m,nTau,numChannel,numThetan,ConstSub);
gt=g(:);
