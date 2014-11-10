function [f,g]=sfun_TensorJ(W,xrfData,xtmData,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0)
global NumSSDlet numChannel NoSelfAbsorption XTMscale gama
global mtol SigMa_XTM SigMa_XRF eX eY 
global LogScale Beta

Beta=1e0;
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,1,NumElement);
MUe_XTM=reshape(MU_e(:,1,1).*gama,1,1,NumElement).*XTMscale;
MU_XTM=sum(reshape(W,m(1),m(2),NumElement).*repmat(MUe,[m(1),m(2),1]),3).*XTMscale;
%%%%% ====================================================================
f=0;
W=reshape(W,mtol,NumElement);
L=reshape(L,length(thetan),nTau+1,mtol);
%%%%% ====================================================================
g=zeros(mtol,NumElement);
for n=1:length(thetan)
    if(LogScale)
        Mt=-log(xtmData(:,n)./I0);
    else
        Mt=xtmData(:,n);
    end
    sum_Tau=0;
    InTens=ones(nTau+1,mtol);
    OutTens=ones(nTau+1,mtol,NumElement);
    for i=1:nTau+1
        XRF_v=zeros(1,numChannel);
        counted_v=[];
        OutTens_d=ones(mtol,NumSSDlet,NumElement);
        index=GlobalInd{n,i};
        if(~isempty(index))
            index_sub=sub2ind(m,index(:,2),index(:,1));
            msub=length(index_sub);
            TempSub=zeros(mtol,NumElement,numChannel);
            for v_count=1:msub
                v=index_sub(v_count);
                v5=SelfInd{n,i,v}{5};
                v7=cell2mat(SelfInd{n,i,v}{7});
                related_v=[v;v5;v7'];
                nonCounted=find(ismember(related_v,counted_v)==0);
                counted_v=[counted_v;related_v(nonCounted)];
                for sub_i=1:length(nonCounted)
                    sub_v=related_v(sub_i);
                    if(isempty(SelfInd{n,i,sub_v}{1}))
                        InTens(i,sub_v)=1;
                    else
                        InTens(i,sub_v)=exp(-sum(sum(W(SelfInd{n,i,sub_v}{1},:).*SelfInd{n,i,sub_v}{3})));
                    end
                    
                    if(~isempty(SelfInd{n,i,sub_v}{2}) & ~NoSelfAbsorption)
                        for d=1:NumSSDlet
                            OutTens_d(sub_v,d,:)=exp(-sum(sum(bsxfun(@times,W(SelfInd{n,i,sub_v}{2}{d},:),reshape(SelfInd{n,i,sub_v}{4}{d},length(SelfInd{n,i,sub_v}{2}{d}),NumElement,NumElement)),1),2));   
                        end
                        OutTens(i,sub_v,:)=sum(OutTens_d(sub_v,:,:),2)/NumSSDlet;
                    end
                end
                v7_tol=[];
                if(~NoSelfAbsorption)
                    temp_sub=zeros(NumElement,numChannel);
                    for i_sub=1:i
                        v7=cell2mat(SelfInd{n,i_sub,v}{7});
                        v7_tol=[v7_tol,v7];
                        OutTens_J=zeros(length(v7),NumElement,NumElement);
                        TempInd=[];
                        for d=1:NumSSDlet
                            TempInd=[length(TempInd)+1:length(SelfInd{n,i_sub,v}{7}{d})+length(TempInd)];
                            if(~isempty(SelfInd{n,i_sub,v}{7}{d}))
                                OutTens_J(TempInd,:,:)=OutTens_J(TempInd,:,:)+bsxfun(@times,OutTens_d(SelfInd{n,i_sub,v}{7}{d},d,:),permute(SelfInd{n,i_sub,v}{8}{d},[3 1 2]));
                            end
                        end
                        OutTens_J=OutTens_J/NumSSDlet;
                        if(~isempty(v7))
                            temp_sub=temp_sub-squeeze(sum(bsxfun(@times,reshape(L(n,i_sub,v7),length(v7),1).*InTens(i_sub,v7)',bsxfun(@times,W(v7,:),OutTens_J)),1))*M;
                        end
                    end
                end
                XRF_v=XRF_v+L(n,i,v)*(InTens(i,v)*reshape(OutTens(i,v,:),1,NumElement).*W(v,:))*M;
                
                if(~isempty(v5) & isempty(v7_tol))
                    TempSub(v,:,:)=-bsxfun(@times,reshape(L(n,i,v5),1,length(v5)).*InTens(i,v5),SelfInd{n,i,v}{6}')*(W(v5,:).*reshape(OutTens(i,v5,:),length(v5),NumElement)*M)...
                        +L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M); % part 3
                elseif(~isempty(v5) & ~isempty(v7_tol))
                    TempSub(v,:,:)=-bsxfun(@times,reshape(L(n,i,v5),1,length(v5)).*InTens(i,v5),SelfInd{n,i,v}{6}')*(W(v5,:).*reshape(OutTens(i,v5,:),length(v5),NumElement)*M)...
                        +temp_sub...
                        +L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M); % part 3
                elseif(isempty(v5) & ~isempty(v7_tol))
                    TempSub(v,:,:)=temp_sub...
                        +L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M); % part 3   
                else
                    TempSub(v,:,:)=L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M); % part 3
                end
            end
            g=g+2*SigMa_XRF((nTau+1)*(n-1)+i)*sum(bsxfun(@times,TempSub,reshape(XRF_v-xrfData{n,i},1,1,numChannel)),3);
            sum_Tau=sum_Tau+SigMa_XRF((nTau+1)*(n-1)+i)*(xrfData{n,i}-XRF_v)*(xrfData{n,i}-XRF_v)';
            Lsub=reshape(L(n,i,:),m(1),m(2));
             if(LogScale)
                count=(nTau+1)*(n-1)+i;
                Rdis=eX'*(MU_XTM.*Lsub)*eY; %% Discrete case
                sum_Tau=sum_Tau+Beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g+reshape(2*Beta*SigMa_XTM(count)*(Rdis-Mt(i)).*repmat(full(Lsub),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]),mtol,NumElement);
            else
                count=(nTau+1)*(n-1)+i;
                Rdis=I0*exp(-eX'*(MU_XTM.*Lsub)*eY);
                sum_Tau=sum_Tau+Beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g-reshape(2*Beta*SigMa_XTM(count)*Rdis*(Rdis-Mt(i)).*repmat(full(Lsub),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]),mtol,NumElement);
            end
        end
    end
    f=f+sum_Tau;
end
g=g(:);
