function [ft,gt]=sfun_full_linear(W,xrfData,M,NumElement,L,GlobalInd,SelfInd,m,nTau)
%%==== XRF objective function in least square form
%%==== Self-absorption is implemented in a tensor-product fashion
%%==== Solve W when the W in exponential term is fixed from the previous iteration
global NumSSDlet numChannel  numThetan % integer: ns nc nt
global NoSelfAbsorption %% bool variable
global W0 InTens OutTens  
mtol=prod(m);
W=repmat(reshape(W,[1 mtol NumElement]),[nTau+1 1 1 numChannel]);
L=reshape(L,numThetan,nTau+1,mtol);
%%%%% ====================================================================
f=zeros(numThetan,1);
g=zeros(numThetan,mtol,NumElement);
for n=1:numThetan
    %%====================================
    temp_v=L(n,:,:).*InTens(n,:,:);
    TempSub=bsxfun(@times,squeeze(bsxfun(@times,temp_v,OutTens(n,:,:,:))),M); % part 3
    XRF_v=sum(sum(TempSub.*W,2),3);
    %%======================================        
    count=(nTau+1)*(n-1)+1:(nTau+1)*n;
    f(n)=sum(sum((squeeze(XRF_v)-squeeze(xrfData(n,:,:))).^2,2),1);
    g(n,:,:)=2*reshape(sum(sum(bsxfun(@times,TempSub, ...
        reshape((squeeze(XRF_v)-squeeze(xrfData(n,:,:))),nTau+1,1,1,numChannel)),1),4),mtol,NumElement);
end
ft=sum(f);
gt=reshape(sum(g,1),mtol*NumElement,1);
