function [ft,gt]=sfun_Tensor_Patrick(W,xrfData,M,NumElement,L,GlobalInd,SelfInd,m,nTau)
%%==== XRF objective function in least square form
%%==== Self-absorption is implemented in a tensor-product fashion
%%==== Solve W when the W in exponential term is fixed from the previous iteration
global NumSSDlet numChannel NoSelfAbsorption numThetan
global SigMa_XRF Wold W0
mtol=prod(m);
W=reshape(W,mtol,NumElement);
Wold=reshape(W0,mtol,NumElement);
L=reshape(L,numThetan,nTau+1,mtol);
%%%%% ====================================================================
f=zeros(numThetan,1);
g=zeros(numThetan,mtol,NumElement);
for n=1:numThetan
    InTens=ones(nTau+1,mtol);
    OutTens=ones(nTau+1,mtol,NumElement);
    OutTens_d=ones(nTau+1,mtol,NumSSDlet,NumElement);
    TempSub=zeros(nTau+1,mtol,NumElement,numChannel);
    XRF_v=zeros(nTau+1,numChannel);
    for i=1:nTau+1
        counted_v=zeros(mtol,1);
        index=GlobalInd{n,i};
        if(~isempty(index))
            index_sub=sub2ind(m,index(:,2),index(:,1));
            msub=length(index_sub);
            for v_count=1:msub
                v=index_sub(v_count);
                if(~isempty(SelfInd{n,i,v}{1}))
                    InTens(i,v)=exp(-sum(sum(Wold(SelfInd{n,i,v}{1},:).*SelfInd{n,i,v}{3})));
                end
                if(~isempty(SelfInd{n,i,v}{2})& ~NoSelfAbsorption)
                    for d=1:NumSSDlet
                        if(~isempty(SelfInd{n,i,v}{2}{d}))
                        OutTens_d(i,v,d,:)=exp(-sum(sum(bsxfun(@times,Wold(SelfInd{n,i,v}{2}{d},:),reshape(SelfInd{n,i,v}{4}{d},length(SelfInd{n,i,v}{2}{d}),NumElement,NumElement)),1),2));
                        end
                    end
                OutTens(i,v,:)=sum(OutTens_d(i,v,:,:),3)/NumSSDlet;
                end
                XRF_v(i,:)=XRF_v(i,:)+L(n,i,v)*(InTens(i,v)*reshape(OutTens(i,v,:),1,NumElement).*W(v,:))*M;
                TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,NumElement,numChannel); % part 3
            end
        end
    end
    count=(nTau+1)*(n-1)+1:(nTau+1)*n;
    f(n)=sum(SigMa_XRF(count).*sum((XRF_v-squeeze(xrfData(n,:,:))).^2,2),1);
    g(n,:,:)=2*reshape(sum(sum(bsxfun(@times,TempSub,repmat(SigMa_XRF(count),[1 1 1 numChannel]).* ...
        reshape((XRF_v-squeeze(xrfData(n,:,:))),nTau+1,1,1,numChannel)),1),4),mtol,NumElement);
end
ft=sum(f);
gt=reshape(sum(g,1),mtol*NumElement,1);
