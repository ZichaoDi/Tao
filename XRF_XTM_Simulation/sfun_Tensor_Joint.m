function [f,g,f_XRF,f_XTM]=sfun_Tensor_Joint(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0)
global NumSSDlet numChannel numThetan NoSelfAbsorption XTMscale
global  SigMa_XTM SigMa_XRF
global LogScale Beta TempBeta 
% load WeightMatrix
f=0;
TempBeta=1;
f_XRF=0;
f_XTM=0;
mtol=prod(m);
W=reshape(W,mtol,NumElement);
L=reshape(full(L),numThetan,nTau+1,m(1),m(2));
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,1,NumElement);
MUe_XTM=reshape(MU_e(:,1,1),1,1,NumElement).*XTMscale;
MU_XTM=sum(reshape(W,m(1),m(2),NumElement).*repmat(MUe,[m(1),m(2),1]),3).*XTMscale;
%%%%% ====================================================================
eX=ones(m(1),1);
eY=ones(m(2),1);
%%%%% ====================================================================
g=zeros(mtol,NumElement);
for n=1:numThetan
    if(LogScale)
        Mt=-log(DisR(:,n)./I0); 
    else
        Mt=DisR(:,n);
    end
    SigMa1=zeros(numChannel,nTau+1);
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
                counted_v=[counted_v;related_v(nonCounted)];
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
                                        reshape(L(sub2ind([numThetan,nTau+1,prod(m)],n,i,sub_v))*InTens(i,sub_v)/NumSSDlet...
                                        *bsxfun(@times,reshape(OutTens_d(i,sub_v,d,:),1,NumElement).*W(sub_v,:),reshape(SelfInd{n,i,v_self}{8}{d}(:,:,sub_v==SelfInd{n,i,v_self}{7}{d}),NumElement,NumElement))*M,[1 1 NumElement numChannel]);
                                end
                            end
                        end
                        OutTens(i,sub_v,:)=sum(OutTens_d(i,sub_v,:,:),3)/NumSSDlet;
                    end
                end
                XRF_v(i,:)=XRF_v(i,:)+L(sub2ind([numThetan,nTau+1,prod(m)],n,i,v))*(InTens(i,v)*reshape(OutTens(i,v,:),1,NumElement).*W(v,:))*M;
                
                if(~isempty(v5))
                    TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(-bsxfun(@times,reshape(L(sub2ind([numThetan,nTau+1,prod(m)],n*ones(size(v5)),i*ones(size(v5)),v5)),1,length(v5)).*InTens(i,v5),SelfInd{n,i,v}{6}')*(W(v5,:).*reshape(OutTens(i,v5,:),length(v5),NumElement)*M)...
                        +L(sub2ind([numThetan,nTau+1,prod(m)],n,i,v))*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,NumElement,numChannel); % part 3
                else
                    TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(full(L(sub2ind([numThetan,nTau+1,prod(m)],n,i,v)))*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,NumElement,numChannel); % part 3
                end
                
            end
            Lsub=full(reshape(L(sub2ind([numThetan,nTau+1,prod(m)],n*ones(1,prod(m)),i*ones(1,prod(m)),1:prod(m))),m(1),m(2)));
            if(LogScale)
                count=(nTau+1)*(n-1)+i;
                Rdis=eX'*(MU_XTM.*Lsub)*eY; %% Discrete case
                f=f+Beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                f_XTM=f_XTM+SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g+reshape(2*Beta*SigMa_XTM(count)*(Rdis-Mt(i)).*repmat(full(Lsub),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]),mtol,NumElement);
            else
                count=(nTau+1)*(n-1)+i;
                Rdis=I0*exp(-eX'*(MU_XTM.*Lsub)*eY);
                f=f+Beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                f_XTM=f_XTM+SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g-reshape(2*Beta*SigMa_XTM(count)*Rdis*(Rdis-Mt(i)).*repmat(full(Lsub),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]),mtol,NumElement);
            end
%              SigMa1(:,i)=diag(diag(-inv(WeightMatrix{n,i}*WeightMatrix{n,i}')))*(XRF_v{n,i}-XRF{n,i})';
        end
    end
    count=(nTau+1)*(n-1)+1:(nTau+1)*n;
    f_XRF=f_XRF+sum(SigMa_XRF(count).*sum((XRF_v-reshape(XRF(n,:,:),nTau+1,numChannel)).^2,2),1);
    f=f+TempBeta*sum(SigMa_XRF(count).*sum((XRF_v-reshape(XRF(n,:,:),nTau+1,numChannel)).^2,2),1);
    g=g+TempBeta*2*reshape(sum(sum(bsxfun(@times,TempSub,repmat(SigMa_XRF(count),[1 1 1 numChannel]).*reshape((XRF_v-reshape(XRF(n,:,:),nTau+1,numChannel)),nTau+1,1,1,numChannel)),1),4),mtol,NumElement);
%     clear InTens OutTens OutTens_d TempSub XRF_v
%     f=f+sum(sum((cat(1,XRF_v{n,:})-cat(1,XRF{n,:})).*SigMa1',2),1);
%     g=g+2*reshape(sum(sum(permute(TempSub,[4 1 2 3]).*repmat(SigMa1,[1,1,mtol,NumElement]),1),2),mtol,NumElement);
% %     for i=1:nTau+1
% %     %     g=g+reshape(reshape(TempSub(i,:,:,:),mtol*NumElement,numChannel)*SigMa1(:,i)+(SigMa(i,:)*reshape(TempSub(i,:,:,:),mtol*NumElement,numChannel)')',mtol,NumElement);
% %         g=g+2*reshape(reshape(TempSub(i,:,:,:),mtol*NumElement,numChannel)*SigMa1(:,i),mtol,NumElement);
% %     
% %     end
end
g=g(:);
