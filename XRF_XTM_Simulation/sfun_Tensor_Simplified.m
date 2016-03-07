function [ft,gt]=sfun_Tensor(W,xrfData,M,NumElement,L,GlobalInd,SelfInd,m,nTau)
%%==== XRF objective function in least square form
%%==== Self-absorption is implemented in a tensor-product fashion
global NumSSDlet numChannel NoSelfAbsorption numThetan
global SigMa_XRF EM MU_e area_xrf

mtol=prod(m);
W=reshape(W,mtol,NumElement);
L=reshape(L,numThetan,nTau+1,mtol);
MU=W*reshape(MU_e(:,1,:),NumElement,NumElement+1);
%%%%% ====================================================================
f=zeros(numThetan,1);
g=zeros(numThetan,mtol,NumElement);
for n=1:numThetan
    InTens=ones(nTau+1,mtol);
    OutTens=ones(nTau+1,mtol,NumElement);
    TempSub=zeros(nTau+1,mtol,NumElement,numChannel);
    XRF_v=zeros(nTau+1,numChannel);
    for i=1:nTau+1
        counted_v=[];
        index=GlobalInd{n,i};
        if(~isempty(index))
            index_sub=sub2ind(m,index(:,2),index(:,1));
            msub=length(index_sub);
            for v_count=1:msub
                v=index_sub(v_count);
                v5=SelfInd{n,i,v}{5};
                if(isempty(SelfInd{n,i,v}{1}))
                    InTens(i,v)=1;
                else
                    InTens(i,v)=exp(-sum(sum(W(SelfInd{n,i,v}{1},:).*SelfInd{n,i,v}{3})));
                end
                
                if(~isempty(SelfInd{n,i,v}{2}) & ~NoSelfAbsorption)
                    OutTens(i,v,:)=exp(-sum(MU(SelfInd{n,i,v}{2},2:NumElement+1),1)*area_xrf(n,i));
                    for d_sub=1:length(SelfInd{n,i,v}{2})
                        v_self=SelfInd{n,i,v}{2}(d_sub);
                        TempSub(i,v_self,:,:)=TempSub(i,v_self,:,:)-...
                            reshape(L(n,i,v)*InTens(i,v)...
                            *repmat(reshape(OutTens(i,v,:),1,NumElement),NumElement,1).*repmat(W(v,:),[NumElement,1]).*squeeze(area_xrf(n,i)*MU_e(:,1,2:NumElement+1))*M,[1 1 NumElement numChannel]);% line 2
                    end
                end
                XRF_v(i,:)=XRF_v(i,:)+L(n,i,v)*(InTens(i,v)*reshape(OutTens(i,v,:),1,NumElement).*W(v,:))*M;
                if(~isempty(v5))
                    TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(-bsxfun(@times,reshape(L(n,i,v5),1,length(v5)).*InTens(i,v5),SelfInd{n,i,v}{6}')*(W(v5,:).*reshape(OutTens(i,v5,:),length(v5),NumElement)*M)... % line 1
                    +L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,NumElement,numChannel); % line 3
                else
                    TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,NumElement,numChannel); % line 3
                end
                
            end
        end
    end
    if(EM)
        XRF_v = XRF_v+1;
        xrfData=xrfData+1;
        f(n)=sum(sum(squeeze(-xrfData(n,:,:)).*log(XRF_v)+XRF_v,1),2);
        g(n,:,:)=1*reshape(sum(sum(bsxfun(@times,TempSub, ...
                reshape((1-squeeze(xrfData(n,:,:))./XRF_v),nTau+1,1,1,numChannel)),1),4),mtol,NumElement);
    else
        f(n)=sum(sum((XRF_v-squeeze(xrfData(n,:,:))).^2,2),1);
        g(n,:,:)=2*reshape(sum(sum(bsxfun(@times,TempSub, ...
                reshape((XRF_v-squeeze(xrfData(n,:,:))),nTau+1,1,1,numChannel)),1),4),mtol,NumElement);
    end
end
ft=sum(f);
gt=reshape(sum(g,1),mtol*NumElement,1);
