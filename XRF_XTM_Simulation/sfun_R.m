function [f, g, shift_y] = sfun_R (W)
%--------------------------------------------------------------
% compute objective function and gradient
%--------------------------------------------------------------
global NumElement M numThetan
global NumSSDlet numChannel NoSelfAbsorption
global SigMa_XRF L GlobalInd SelfInd m nTau XRF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=0;
mtol=prod(m);
W=reshape(W,mtol,NumElement);
shift_y=zeros(numThetan,nTau+1,numChannel);
%%%%% ====================================================================
g=zeros(mtol,NumElement);
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
                                        *bsxfun(@times,reshape(OutTens_d(i,sub_v,d,:),1,NumElement).*W(sub_v,:),reshape(SelfInd{n,i,v_self}{8}{d}(:,:,sub_v==SelfInd{n,i,v_self}{7}{d}),NumElement,NumElement))*M,[1 1 NumElement numChannel]);
                                end
                            end
                        end
                        OutTens(i,sub_v,:)=sum(OutTens_d(i,sub_v,:,:),3)/NumSSDlet;
                    end
                end
                XRF_v(i,:)=XRF_v(i,:)+L(n,i,v)*(InTens(i,v)*reshape(OutTens(i,v,:),1,NumElement).*W(v,:))*M;
                
                if(~isempty(v5))
                    TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(-bsxfun(@times,reshape(L(n,i,v5),1,length(v5)).*InTens(i,v5),SelfInd{n,i,v}{6}')*(W(v5,:).*reshape(OutTens(i,v5,:),length(v5),NumElement)*M)...
                        +L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,NumElement,numChannel); % part 3
                else
                    TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,NumElement,numChannel); % part 3
                end
                
            end
        end
    end
    count=(nTau+1)*(n-1)+1:(nTau+1)*n;
    shift_y(n,:,:)=XRF_v;
    f=f+sum(SigMa_XRF(count).*sum((XRF_v-squeeze(XRF(n,:,:))).^2,2),1);
    g=g+2*reshape(sum(sum(bsxfun(@times,TempSub,repmat(SigMa_XRF(count),[1 1 1 numChannel]).*reshape((XRF_v-squeeze(XRF(n,:,:))),nTau+1,1,1,numChannel)),1),4),mtol,NumElement);
end
g=g(:);

