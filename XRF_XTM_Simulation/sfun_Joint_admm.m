function [f,g,f_XRF,f_XTM]=sfun_Joint_admm(W,xrfData,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,m,nTau)
%%==== XRF objective function in least square form
%%==== Self-absorption is implemented in a tensor-product fashion
%%==== Solve W when the W in exponential term is fixed from the previous iteration
global NumSSDlet numChannel numThetan I0 
global  SigMa_XTM SigMa_XRF W0 
global InTens OutTens MU_XTM
global LogScale Beta TempBeta NoSelfAbsorption XTMscale
mtol=prod(m);
W=reshape(W,mtol,NumElement);
L=reshape(L,numThetan,nTau+1,mtol);
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,NumElement);
MUe_XTM=MUe*XTMscale;
MU_XTM=sum(W.*repmat(MUe,[mtol,1]),2)'*XTMscale;
%%%%% ====================================================================
f_XRF=0;
f_XTM=0;
g=zeros(numThetan,mtol,NumElement);
for n=1:numThetan
    if(LogScale)
        Mt=-log(DisR(:,n)./I0);
    else
        Mt=DisR(:,n);
    end
    %%====================================
    temp_v=L(n,:,:).*InTens(n,:,:);
    TempSub=repmat(squeeze(temp_v),[1,1,NumElement numChannel]).*repmat(squeeze(OutTens(n,:,:,:)),[1,1,1,numChannel]).*repmat(reshape(M,[1,1 NumElement numChannel]),[nTau+1,mtol,1,1]); % part 3
    XRF_v=sum(sum(squeeze(TempSub).*repmat(reshape(W,[1 mtol NumElement]),[nTau+1 1 1 numChannel]),2),3);
    %%======================================        
    Lsub=squeeze(L(n,:,:));
    if(LogScale)
        Rdis=sum(repmat(MU_XTM,[nTau+1,1]).*Lsub,2); %% Discrete case
        f_XTM=f_XTM+norm(Rdis-Mt)^2;
        g_XTM(n,:,:)=reshape(sum(2*repmat(Rdis-Mt,[1,mtol,NumElement]).*repmat(full(Lsub),[1,1,NumElement]).*repmat(reshape(MUe_XTM,[1,1,NumElement]),[nTau+1,mtol,1]),1),1,mtol,NumElement);
    else
        Rdis=I0*exp(-sum(sum(repmat(reshape(MU_XTM,[1,m(1),m(2)]),[nTau+1,1,1]).*Lsub,2),3)); %% Discrete case
        f_XTM=f_XTM+norm(Rdis-Mt)^2;
    end
    count=(nTau+1)*(n-1)+1:(nTau+1)*n;
    f_XRF=f_XRF+sum(SigMa_XRF(count).*sum((squeeze(XRF_v)-squeeze(xrfData(n,:,:))).^2,2),1);
    g_XRF(n,:,:)=2*reshape(sum(sum(bsxfun(@times,TempSub,repmat(SigMa_XRF(count),[1 1 1 numChannel]).* ...
        reshape((squeeze(XRF_v)-squeeze(xrfData(n,:,:))),nTau+1,1,1,numChannel)),1),4),1,mtol,NumElement);

end
f=TempBeta*f_XRF+Beta*f_XTM;
g=reshape(sum(TempBeta*g_XRF+Beta*g_XTM,1),mtol*NumElement,1);
