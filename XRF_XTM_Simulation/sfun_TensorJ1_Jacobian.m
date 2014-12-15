function [f,g]=sfun_TensorJ1(W,xrfData,xtmData,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0)
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
Jacob1=[];
Jacob2=[];
for n=1:length(thetan)
    if(LogScale)
        Mt=-log(xtmData(:,n)./I0);
    else
        Mt=xtmData(:,n);
    end
    sum_Tau=0;
    InTens=ones(nTau+1,mtol);
    OutTens=ones(nTau+1,mtol,NumElement);
    OutTens_d=ones(nTau+1,mtol,NumSSDlet,NumElement);
    for i=1:nTau+1
        JacobSub1=zeros(mtol,NumElement,numChannel);
        XRF_v=zeros(1,numChannel);
        counted_v=[];
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
                            OutTens_d(i,sub_v,d,:)=exp(-sum(sum(bsxfun(@times,W(SelfInd{n,i,sub_v}{2}{d},:),reshape(SelfInd{n,i,sub_v}{4}{d},length(SelfInd{n,i,sub_v}{2}{d}),NumElement,NumElement)),1),2));
                        end
                        OutTens(i,sub_v,:)=sum(sum(OutTens_d(i,sub_v,:,:),1),3)/NumSSDlet;
                    end
                end
                temp_sub=zeros(NumElement,numChannel);
                if(~NoSelfAbsorption)
                    for i_sub=1:i
                        v7=cell2mat(SelfInd{n,i_sub,v}{7});
                        if(~isempty(v7))
                        OutTens_J=sum(bsxfun(@times,permute(OutTens_d(i_sub,v7,:,:),[2 1 4 3]),permute(cat(3,SelfInd{n,i_sub,v}{8}{:}),[3 1 2])),4)/(NumSSDlet^2);
                        temp_sub=temp_sub-squeeze(sum(bsxfun(@times,reshape(L(n,i_sub,v7),length(v7),1).*InTens(i_sub,v7)',bsxfun(@times,W(v7,:),OutTens_J)),1))*M;
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Temp1=cell(i,1);
%                     Temp1=SelfInd{n,1:i,v}
%                     Temp=cellfun(@(x)x{7},SelfInd{n,1:i,v},'UniformOutput',false)
%                     pause;
%                     
%  gLinear=cell2mat(Temp1);
%  if(~isempty(gLinear))
%      maxLength=max(cellfun(@(x)numel(x),gLinear));
%      out=cell2mat(cellfun(@(x)cat(2,x,zeros(1,maxLength-length(x))),gLinear,'UniformOutput',false))
% pause;
%                     TempInd=horzcat(SelfInd{n,1:i,v});
%                     v7=cat(2,TempInd{7,:});
%                     if(~isempty(v7))
%                         TempInd=cat(3,horzcat(TempInd{8,:}));
%                         OutTens_J=sum(bsxfun(@times,reshape(permute(OutTens_d(1:i,v7,:,:),[1 3 2]),length(v7),1,NumElement,NumSSDlet),permute(TempInd,[3 1 2])),4)/(NumSSDlet^2);
%                         temp_sub=temp_sub-squeeze(sum(bsxfun(@times,reshape(L(n,i_sub,v7),length(v7),1).*InTens(i_sub,v7)',bsxfun(@times,W(v7,:),OutTens_J)),1))*M;
%                     end
%                     pause;
%  end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                XRF_v=XRF_v+L(n,i,v)*(InTens(i,v)*reshape(OutTens(i,v,:),1,NumElement).*W(v,:))*M;
                
                if(~isempty(v5))
                    TempSub(v,:,:)=-bsxfun(@times,reshape(L(n,i,v5),1,length(v5)).*InTens(i,v5),SelfInd{n,i,v}{6}')*(W(v5,:).*reshape(OutTens(i,v5,:),length(v5),NumElement)*M)+ ...
                        temp_sub+L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M); % part 3
                else
                    TempSub(v,:,:)=temp_sub+L(n,i,v)*InTens(i,v)*bsxfun(@times,squeeze(OutTens(i,v,:)),M); % part 3
                end
            JacobSub1(v,:,:)=JacobSub1(v,:,:)+TempSub(v,:,:);
            end
            g=g+2*SigMa_XRF((nTau+1)*(n-1)+i)*sum(bsxfun(@times,TempSub,reshape(XRF_v-xrfData{n,i},1,1,numChannel)),3);
            sum_Tau=sum_Tau+SigMa_XRF((nTau+1)*(n-1)+i)*(xrfData{n,i}-XRF_v)*(xrfData{n,i}-XRF_v)';
            Lsub=reshape(L(n,i,:),m(1),m(2));
            if(LogScale)
                count=(nTau+1)*(n-1)+i;
                Rdis=eX'*(MU_XTM.*Lsub)*eY; %% Discrete case
                sum_Tau=sum_Tau+Beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g+reshape(2*Beta*SigMa_XTM(count)*(Rdis-Mt(i)).*repmat(full(Lsub),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]),mtol,NumElement);
            JacobSub2=2*Beta*SigMa_XTM(count)*repmat(full(Lsub),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]);
            else
                count=(nTau+1)*(n-1)+i;
                Rdis=I0*exp(-eX'*(MU_XTM.*Lsub)*eY);
                sum_Tau=sum_Tau+Beta*SigMa_XTM(count)*(Rdis-Mt(i))^2;
                g=g-reshape(2*Beta*SigMa_XTM(count)*Rdis*(Rdis-Mt(i)).*repmat(full(Lsub),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]),mtol,NumElement);
             JacobSub2=-2*Beta*SigMa_XTM(count)*Rdis*repmat(full(Lsub),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]);
            end
             Jacob1=[Jacob1;reshape(JacobSub1,m(1)*m(2)*NumElement,numChannel)'];
            Jacob2=[Jacob2;JacobSub2(:)'];
        end
    end
    f=f+sum_Tau;
end
%%=============================================== Measure ill-conditioness
% figure,subplot(2,1,1),spy(Jacob1), subplot(2,1,2),spy(Jacob2);
size(Jacob1)
size(Jacob2)
tol=eps;
J1=rank(Jacob1,tol)
% k1=cond(Jacob1);
% [U1,S1,V1]=svd(Jacob1);
J2=rank(Jacob2,tol)
% k2=cond(Jacob2);
% [U2,S2,V2]=svd(Jacob2);
% [U,S,V]=svd([Jacob1;Jacob2]);
J3=rank([Jacob1;Jacob2],tol)
% fprintf('Cond1=%d, Rank1=%d, Cond2=%d, Rank2=%d, cond=%d, minS=%d, rank=%d \n',k1,J1,k2,J2,cond([Jacob1;Jacob2]),min(diag(S)),J3);
%%=====================================================================
g=g(:);
