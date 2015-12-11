function [f,g]=sfun_linearized(DW,xrfData,M,NumElement,L,m,nTau)
%%==== XRF objective function in least square form
%%==== Aggregate Self-absorption, beam attanuation and W as a single variable and minize for it.
global numChannel  numThetan
mtol=prod(m);
L=reshape(L,numThetan,nTau+1,mtol);
DW=reshape(DW,numThetan,nTau+1,mtol,NumElement);
%%%%% ====================================================================
XRF_v=squeeze(sum(sum(bsxfun(@times,bsxfun(@times,L,DW),M),3),4));
TempSub=bsxfun(@times,L,M); 
%%=================== forwad part only
% f=sum(sum(XRF_v,1),2);
% g=reshape(TempSub,numThetan*(nTau+1)*mtol*NumElement,numChannel);
% return;
%%==============================
f=(XRF_v-xrfData).^2;
f=sum(f(:));
g=2*sum(bsxfun(@times,TempSub, ...
        reshape(XRF_v-xrfData,[numThetan,nTau+1,1,1,numChannel])),5);
g=g(:);
