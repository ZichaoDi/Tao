function [f,g,f_XRF,f_XTM]=sfun_Joint_admm(W,xrfData,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,m,nTau)
%%==== XRF objective function in least square form
%%==== Self-absorption is implemented in a tensor-product fashion
%%==== Solve W when the W in exponential term is fixed from the previous iteration
global NumSSDlet numChannel numThetan I0 W0 
global InTens OutTens MU_XTM
global LogScale Beta TempBeta NoSelfAbsorption
mtol=prod(m);
W=repmat(reshape(W,[1 mtol NumElement]),[nTau+1 1 1 numChannel]);
L=reshape(L,numThetan,nTau+1,mtol);
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,NumElement);
MUe_XTM=MUe*1;
MU_XTM=sum(bsxfun(@times,squeeze(W(1,:,:,1)),MUe),2)';
%%%%% ====================================================================
f_XRF=zeros(numThetan,1);
f_XTM=zeros(numThetan,1);
g=zeros(numThetan,mtol,NumElement);
if(LogScale)
    Mt=-log(DisR./I0);
else
    Mt=DisR;
end
for n=1:numThetan
    %%====================================
    temp_v=L(n,:,:).*InTens(n,:,:);
    TempSub=bsxfun(@times,squeeze(bsxfun(@times,temp_v,OutTens(n,:,:,:))),M); % part 3
    XRF_v=sum(sum(TempSub.*W,2),3);
    %%======================================        
    Lsub=squeeze(L(n,:,:));
    if(LogScale)
        Rdis=sum(repmat(MU_XTM,[nTau+1,1]).*Lsub,2); %% Discrete case
        f_XTM(n)=norm(Rdis-Mt(:,n))^2;
        g_XTM(n,:,:)=reshape(sum(2*bsxfun(@times,bsxfun(@times,Rdis-Mt(:,n),Lsub),reshape(MUe_XTM,[1,1,NumElement])),1),1,mtol,NumElement);
    end
    count=(nTau+1)*(n-1)+1:(nTau+1)*n;
    f_XRF(n)=sum(sum((squeeze(XRF_v)-squeeze(xrfData(n,:,:))).^2,2),1);
    g_XRF(n,:,:)=2*reshape(sum(sum(bsxfun(@times,TempSub, ...
        reshape((squeeze(XRF_v)-squeeze(xrfData(n,:,:))),nTau+1,1,1,numChannel)),1),4),1,mtol,NumElement);

end
f_XRF=sum(f_XRF);
f_XTM=sum(f_XTM);
f=TempBeta*f_XRF+Beta*f_XTM;
g=reshape(sum(TempBeta*g_XRF+Beta*g_XTM,1),mtol*NumElement,1);
