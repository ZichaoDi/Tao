function [ft,gt,Jacob1,Jacob2,Jacob3]=sfun_Tensor(W,xrfData,M,NumElement,L,GlobalInd,SelfInd,m,nTau)
%%==== XRF objective function in least square form
%%==== Self-absorption is implemented in a tensor-product fashion
global NumSSDlet numChannel NoSelfAbsorption numThetan
global SigMa_XRF EM

mtol=prod(m);
W=reshape(W,mtol,NumElement);
L=reshape(L,numThetan,nTau+1,mtol);
%%%%% ====================================================================
f=zeros(numThetan,1);
g=zeros(numThetan,mtol,NumElement);
Jacob1=zeros(numThetan,nTau+1,mtol,NumElement,numChannel);
Jacob2=zeros(numThetan,nTau+1,mtol,NumElement,numChannel);
Jacob3=zeros(numThetan,nTau+1,mtol,NumElement,numChannel);
for n=1:numThetan
    InTens=ones(nTau+1,mtol);
    OutTens=ones(nTau+1,mtol,NumElement);
    OutTens_d=ones(nTau+1,mtol,NumSSDlet,NumElement);
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
                    
                    if(~isempty(SelfInd{n,i,sub_v}{2}) & ~NoSelfAbsorption)
                        for d=1:NumSSDlet
                            if(~isempty(SelfInd{n,i,sub_v}{2}{d}))
                                OutTens_d(i,sub_v,d,:)=exp(-sum(sum(bsxfun(@times,W(SelfInd{n,i,sub_v}{2}{d},:),reshape(SelfInd{n,i,sub_v}{4}{d},length(SelfInd{n,i,sub_v}{2}{d}),NumElement,NumElement)),1),2));
                                for d_sub=1:length(SelfInd{n,i,sub_v}{2}{d})
                                    v_self=SelfInd{n,i,sub_v}{2}{d}(d_sub);
                                    TempSub(i,v_self,:,:)=TempSub(i,v_self,:,:)-...
                                        reshape(L(n,i,sub_v)*InTens(i,sub_v)/NumSSDlet...
                                        *bsxfun(@times,reshape(OutTens_d(i,sub_v,d,:),1,NumElement).*W(sub_v,:),reshape(SelfInd{n,i,v_self}{8}{d}(:,:,sub_v==SelfInd{n,i,v_self}{7}{d}),NumElement,NumElement))*M,[1 1 NumElement numChannel]);% line 2
                                    Jacob1(n,i,v_self,:,:)=Jacob1(n,i,v_self,:,:)-...
                                        reshape(L(n,i,sub_v)*InTens(i,sub_v)/NumSSDlet...
                                        *bsxfun(@times,reshape(OutTens_d(i,sub_v,d,:),1,NumElement).*W(sub_v,:),reshape(SelfInd{n,i,v_self}{8}{d}(:,:,sub_v==SelfInd{n,i,v_self}{7}{d}),NumElement,NumElement))*M,[1,1 1 NumElement numChannel]);
                                end
                            end
                        end
                        OutTens(i,sub_v,:)=sum(OutTens_d(i,sub_v,:,:),3)/NumSSDlet;
                    end
                end
                XRF_v(i,:)=XRF_v(i,:)+L(n,i,v)*(InTens(i,v)*reshape(OutTens(i,v,:),1,NumElement).*W(v,:))*M;
                if(~isempty(v5))
                    TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(-bsxfun(@times,reshape(L(n,i,v5),1,length(v5)).*InTens(i,v5),SelfInd{n,i,v}{6}')*(W(v5,:).*reshape(OutTens(i,v5,:),length(v5),NumElement)*M)... % line 1
                    +L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,NumElement,numChannel); % line 3
                    Jacob2(n,i,v,:,:)=Jacob2(n,i,v,:,:)+1*reshape(-bsxfun(@times,reshape(L(n,i,v5),1,length(v5)).*InTens(i,v5),SelfInd{n,i,v}{6}')*(W(v5,:).*reshape(OutTens(i,v5,:),length(v5),NumElement)*M),1,1,1,NumElement,numChannel); % line 3
                    Jacob3(n,i,v,:,:)=Jacob3(n,i,v,:,:)+1*reshape(L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,1,NumElement,numChannel);

                else
                    TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,NumElement,numChannel); % line 3
                    Jacob3(n,i,v,:,:)=Jacob3(n,i,v,:,:)+1*reshape(L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,1,NumElement,numChannel); 
                end
                
            end
        end
    end
    if(EM)
        XRF_v = XRF_v+1;
        xrfData=xrfData+1;
        f(n)=sum(sum(squeeze(xrfData(n,:,:)).*log(XRF_v)-XRF_v,1),2);
        g(n,:,:)=2*reshape(sum(sum(bsxfun(@times,TempSub, ...
                reshape((1-squeeze(xrfData(n,:,:))./XRF_v),nTau+1,1,1,numChannel)),1),4),mtol,NumElement);
    else
        f(n)=sum(sum((XRF_v-squeeze(xrfData(n,:,:))).^2,2),1);
        g(n,:,:)=2*reshape(sum(sum(bsxfun(@times,TempSub, ...
                reshape((XRF_v-squeeze(xrfData(n,:,:))),nTau+1,1,1,numChannel)),1),4),mtol,NumElement);
    end
end
ft=sum(f);
gt=reshape(sum(g,1),mtol*NumElement,1);
