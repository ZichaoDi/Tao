function [f,g]=sfun_Tensor_Joint(W,xrfData,xtmData,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0)
global NumSSDlet numChannel NoSelfAbsorption XTMscale gama
global mtol SigMa_XTM SigMa_XRF eX eY
global LogScale Beta EmptyBeam
%load WeightMatrix
f=0;
W=reshape(W,mtol,NumElement);
L=reshape(L,length(thetan),nTau+1,mtol);
Beta=1e0;
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,1,NumElement);
MUe_XTM=reshape(MU_e(:,1,1).*gama,1,1,NumElement).*XTMscale;
MU_XTM=sum(reshape(W,m(1),m(2),NumElement).*repmat(MUe,[m(1),m(2),1]),3).*XTMscale;
%%%%% ====================================================================
%%%%% ====================================================================
g=zeros(mtol,NumElement);
for n=1:length(thetan)
    if(LogScale)
        Mt=-log(xtmData(:,n)./I0);
    else
        Mt=xtmData(:,n);
    end
    SigMa1=zeros(numChannel,nTau+1);
    SigMa=SigMa1';
    InTens=ones(nTau+1,mtol);
    OutTens=ones(nTau+1,mtol,NumElement);
    OutTens_d=ones(nTau+1,mtol,NumSSDlet,NumElement);
    TempSub=zeros(nTau+1,mtol,NumElement,numChannel);
    for i=1:nTau+1
        XRF_v{n,i}=zeros(1,numChannel);
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
                            end
                        end
                        OutTens(i,sub_v,:)=sum(OutTens_d(i,sub_v,:,:),3)/NumSSDlet;
                    end
                end
                XRF_v{n,i}=XRF_v{n,i}+L(n,i,v)*(InTens(i,v)*reshape(OutTens(i,v,:),1,NumElement).*W(v,:))*M;
                if(~isempty(v5))
                    TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(-bsxfun(@times,reshape(L(n,i,v5),1,length(v5)).*InTens(i,v5),SelfInd{n,i,v}{6}')*(W(v5,:).*reshape(OutTens(i,v5,:),length(v5),NumElement)*M)...
                        +L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,NumElement,numChannel); % part 3
                else
                    TempSub(i,v,:,:)=TempSub(i,v,:,:)+1*reshape(L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M),1,1,NumElement,numChannel); % part 3
                end
                
            end
            Lsub=reshape(L(n,i,:),m(1),m(2));
            if(LogScale)
                count=(nTau+1)*(n-1)+i;
                Rdis=eX'*(MU_XTM.*Lsub)*eY; %% Discrete case
                f=f+Beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g+reshape(2*Beta*SigMa_XTM(count)*(Rdis-Mt(i)).*repmat(full(Lsub),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]),mtol,NumElement);
            else
                count=(nTau+1)*(n-1)+i;
                Rdis=I0*exp(-eX'*(MU_XTM.*Lsub)*eY);
                f=f+Beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g-reshape(2*Beta*SigMa_XTM(count)*Rdis*(Rdis-Mt(i)).*repmat(full(Lsub),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]),mtol,NumElement);
            end
         % SigMa1(:,i)=eye(size(diag(diag(-inv(WeightMatrix{n,i}*WeightMatrix{n,i}')))))*(XRF_v{n,i}-xrfData{n,i})';
        end
    end
    if(~NoSelfAbsorption)
        UsedBeam=setdiff([1:nTau+1],EmptyBeam(2,EmptyBeam(1,:)==n));
        for i_sub_count=1:length(UsedBeam)
            i_sub=UsedBeam(i_sub_count);
            for v=1:mtol
                v7=cell2mat(SelfInd{n,i_sub,v}{7});
                if(~isempty(v7))
                    OutTens_J=zeros(length(v7),NumElement,NumElement);
                    TempInd=0;
                    for d=1:NumSSDlet
                        if(~isempty(SelfInd{n,i_sub,v}{7}{d}))
                            TempInd=[TempInd(end)+1:length(SelfInd{n,i_sub,v}{7}{d})+TempInd(end)];
                            [~,mx,my,mz]=size(OutTens_d(i_sub,SelfInd{n,i_sub,v}{7}{d},d,:));
                            OutTens_J(TempInd,:,:)=OutTens_J(TempInd,:,:)+bsxfun(@times,reshape(OutTens_d(i_sub,SelfInd{n,i_sub,v}{7}{d},d,:)...
                                ,mx,my,mz),permute(reshape(SelfInd{n,i_sub,v}{8}{d},NumElement,NumElement,length(SelfInd{n,i_sub,v}{7}{d})),[3 1 2]));
                        end
                    end
                    OutTens_J=OutTens_J/NumSSDlet;
                    TempSub(i_sub,v,:,:)=TempSub(i_sub,v,:,:)-reshape(squeeze(sum(bsxfun(@times,reshape(L(n,i_sub,v7),length(v7),1).*InTens(i_sub,v7)'...
                        ,bsxfun(@times,W(v7,:),permute(OutTens_J,[1 3 2]))),1))'*M,1,1,NumElement,numChannel);
                end
            end
        end
    end
    count=(nTau+1)*(n-1)+1:(nTau+1)*n;
    f=f+sum(SigMa_XRF(count).*sum((cat(1,XRF_v{n,:})-cat(1,xrfData{n,:})).^2,2),1);
g=g+2*reshape(sum(sum(bsxfun(@times,TempSub,repmat(SigMa_XRF(count),[1 1 1 numChannel]).*reshape((cat(1,XRF_v{n,:})-cat(1,xrfData{n,:})),nTau+1,1,1,numChannel)),1),4),mtol,NumElement);

% f=f+sum(sum((cat(1,XRF_v{n,:})-cat(1,xrfData{n,:})).*SigMa1',2),1);
% g=g+2*reshape(sum(sum(permute(TempSub,[4 1 2 3]).*repmat(SigMa1,[1,1,mtol,NumElement]),1),2),mtol,NumElement);
% for i=1:nTau+1
% %     g=g+reshape(reshape(TempSub(i,:,:,:),mtol*NumElement,numChannel)*SigMa1(:,i)+(SigMa(i,:)*reshape(TempSub(i,:,:,:),mtol*NumElement,numChannel)')',mtol,NumElement);
%     g=g+2*reshape(reshape(TempSub(i,:,:,:),mtol*NumElement,numChannel)*SigMa1(:,i),mtol,NumElement);
% 
% end
end
g=g(:);