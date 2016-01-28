function [ft,gt]=sfun_half_linear(W,xrfData,M,NumElement,L,GlobalInd,SelfInd,m,nTau)
%%==== XRF objective function in least square form
%%==== Self-absorption is implemented in a tensor-product fashion
global NumSSDlet numChannel NoSelfAbsorption numThetan
global SigMa_XRF OutTens

mtol=prod(m);
W=reshape(W,mtol,NumElement);
L=reshape(L,numThetan,nTau+1,mtol);
%%%%% ====================================================================
f=zeros(numThetan,1);
g=zeros(numThetan,mtol,NumElement);
for n=1:numThetan
    InTens=ones(nTau+1,mtol);
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
                v5=SelfInd{n,i,v}{5};
                v7=cell2mat(SelfInd{n,i,v}{7});
                related_v=unique([v;v5;v7']);
                nonCounted=find(ismember(related_v,counted_v)==0);
                counted_v(related_v(nonCounted))=related_v(nonCounted);
                for sub_i=1:length(nonCounted)
                    sub_v=related_v(nonCounted);
                    sub_v=sub_v(sub_i);
                    if(isempty(SelfInd{n,i,sub_v}{1}))
                        InTens(i,sub_v)=1;
                    else
                        InTens(i,sub_v)=exp(-sum(sum(W(SelfInd{n,i,sub_v}{1},:).*SelfInd{n,i,sub_v}{3})));
                    end
                    
                end
                XRF_v(i,:)=XRF_v(i,:)+L(n,i,v)*(InTens(i,v)*reshape(OutTens(n,i,v,:),1,NumElement).*W(v,:))*M;
                
                if(~isempty(v5))
                    TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(-bsxfun(@times,reshape(L(n,i,v5),1,length(v5)).*InTens(i,v5),SelfInd{n,i,v}{6}')*(W(v5,:).*reshape(OutTens(n,i,v5,:),length(v5),NumElement)*M)...
                        +L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(n,i,v,:)),M),1,1,NumElement,numChannel); % part 3
                else
                    TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(n,i,v,:)),M),1,1,NumElement,numChannel); % part 3
                end
                
            end
        end
    end
    count=(nTau+1)*(n-1)+1:(nTau+1)*n;
    f(n)=sum(SigMa_XRF(count)'.*sum((XRF_v-squeeze(xrfData(n,:,:))).^2,2),1);
    g(n,:,:)=2*reshape(sum(sum(bsxfun(@times,TempSub,repmat(SigMa_XRF(count)',[1 1 1 numChannel]).* ...
        reshape((XRF_v-squeeze(xrfData(n,:,:))),nTau+1,1,1,numChannel)),1),4),mtol,NumElement);
end
ft=sum(f);
gt=reshape(sum(g,1),mtol*NumElement,1);
