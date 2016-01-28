function [f,g]=sfun_linearized_single(DW,xrfData,M,NumElement,L,m)
%%==== XRF objective function in least square form
%%==== Aggregate Self-absorption, beam attanuation and W as a single variable and minimize for it.
global numChannel  numThetan Tik penalty
mtol=prod(m);
L=reshape(L,mtol,1);
DW=reshape(DW,mtol,NumElement);
%%%%% ==========================
XRF_v=squeeze(M(1,:,:))'*DW'*L;
TempSub=bsxfun(@times,L,M(1,:,:)); 
%%==============================
f=(XRF_v-xrfData).^2;
f=sum(f(:));
g=2*sum(bsxfun(@times,TempSub, ...
        reshape(XRF_v-xrfData,[1,numChannel])),2);
g=g(:);
if(penalty)
    W=DW(:);
    lambda=1e-9;
    f=f+lambda*(norm(Tik*W))^2;
    Wdiff=W(1:end-1)-W(2:end);
    g=g+lambda.*(-2.*[0;Wdiff]+2.*[Wdiff;0]);
end
