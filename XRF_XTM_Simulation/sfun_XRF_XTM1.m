function [f,g]=sfun_XRF_XTM(W,xrfData,xtmData,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,LocalInd,L_after,thetan,m,nTau,I0)
global BeforeEmit NumSSDlet NoSelfAbsorption testind XTMscale gama
global SigMa_XTM SigMa_XRF
global W0 VarInd LogScale
f=0;
Wtol=W0;Wtol(VarInd)=W;
W=reshape(Wtol,m(1),m(2),NumElement);
beta=1e0;
alpha=beta;
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,1,NumElement);
MU=sum(W.*repmat(MUe,[m(1),m(2),1]),3);
MUe_XTM=reshape(MU_e(:,1,1).*gama,1,1,NumElement).*XTMscale;
MU_XTM=sum(W.*repmat(MUe,[m(1),m(2),1]),3).*XTMscale;
%%%%% =================== Attenuation Matrix at flourescence energy (Corrected Attenuation)
MU_after=cell(NumElement,1);
for i=1:NumElement
    MU_after{i}=sum(W.*repmat(reshape(MU_e(:,1,i+1),1,1,NumElement),[m(1),m(2),1]),3);
end
%%%%% ====================================================================
e1=ones(m(1),1);
g=zeros(m(1),m(2),NumElement);
Jacob1=[];
Jacob2=[];
SelfInd1=cell(m(1),m(2));
for im=1:m(1)
    for jm=1:m(2)
        SelfInd1{im,jm}=cell(nTau+1,1);
        for i=1:nTau+1
            SelfInd1{im,jm}{i}=cell(5,1);
        end
    end
end
clear i j t
ttt=0;
for n=1:length(thetan)
    SelfInd=SelfInd1;
    if(LogScale)
        Mt=-log(xtmData(:,n)./I0);
    else
        Mt=xtmData(:,n);
    end
    sum_Tau=0;
    for i=1:nTau+1
        JacobSub1=zeros(m(1),m(2),NumElement,numChannel);
        index=GlobalInd{n,i};
        if(~isempty(index))
            ttt=ttt+1;
            L=Ltol{n,i};
            RM=cell(m(1),m(2));
            xrfSub=zeros(1,numChannel);
            TempSub=zeros(NumElement,size(index,1),numChannel);
            I_incident=zeros(size(index,1),1);
            temp_d=zeros(NumElement,m(1),m(2),NumElement);
            for j=1:size(index,1)
                if(j==1)
                    I_incident(j)=1;
                    temp_sum=0;
                else
                    temp_sum=temp_sum+L(index(j-1,2),index(j-1,1))*MU(index(j-1,2),index(j-1,1));
                    I_incident(j)=exp(-temp_sum);
                end
                Wsub=reshape(W(index(j,2),index(j,1),:),[NumElement,1]);
                %% Self-absorption
                
                BeforeEmit=0;
                I_after=ones(NumElement,1);
                Ind_after=[];
                if(NoSelfAbsorption)
                    NumSSDlet=1;%% Turn off self-absorption
                else
                    I_after=0*I_after;
                    
                    for SSDi=1:NumSSDlet
                        index_after=LocalInd{n,i,index(j,2),index(j,1),SSDi};
                        Lvec_after=L_after{n,i,index(j,2),index(j,1),SSDi};
                        Ind_after=unique([Ind_after;index_after],'rows');
                        LinearInd=sub2ind([m(1),m(2)],index_after(:,2),index_after(:,1));
                        for tsub=1:NumElement
                            temp_after=sum(Lvec_after.*MU_after{tsub}(LinearInd)); %% Attenuation of Flourescent energy emitted from current pixel
                            for si=1:size(index_after,1)
                                temp_d(:,index_after(si,2),index_after(si,1),tsub)=temp_d(:,index_after(si,2),index_after(si,1),tsub)+exp(-temp_after)*reshape(MU_e(:,1,tsub+1),NumElement,1).*Lvec_after(si);
                                if(ismac)
                                    TempIs=ismember(index(j,end:-1:1),SelfInd{index_after(si,2),index_after(si,1)}{i}{1},'rows');
                                else
                                    TempIs=ismember(index(j,end:-1:1),SelfInd{index_after(si,2),index_after(si,1)}{i}{1},'rows','legacy');
                                end
                                if(~TempIs)                                  
                                    SelfInd{index_after(si,2),index_after(si,1)}{i}{1}=[SelfInd{index_after(si,2),index_after(si,1)}{i}{1};index(j,end:-1:1)];%% assign downstream index to pixel(oppsite)
                                    SelfInd{index_after(si,2),index_after(si,1)}{i}{2}=[SelfInd{index_after(si,2),index_after(si,1)}{i}{2};L(index(j,2),index(j,1))];
                                    SelfInd{index_after(si,2),index_after(si,1)}{i}{3}=[SelfInd{index_after(si,2),index_after(si,1)}{i}{3};reshape(W(index(j,2),index(j,1),:),1,NumElement)];
                                    SelfInd{index_after(si,2),index_after(si,1)}{i}{4}=[SelfInd{index_after(si,2),index_after(si,1)}{i}{4};I_incident(j)];
                                end
                            end
                            I_after(tsub)=I_after(tsub)+exp(-temp_after)/NumSSDlet;
                        end %% End loop for existing fluorescence energy from current pixel
                    end %% End loop for each SSD detector let
                    for tt=1:size(Ind_after,1)
                        if(isempty(SelfInd{Ind_after(tt,2),Ind_after(tt,1)}{i}{5}))
                            SelfInd{Ind_after(tt,2),Ind_after(tt,1)}{i}{5}=cell(NumElement,1);
                        end
                        for tsub=1:NumElement
                            SelfInd{Ind_after(tt,2),Ind_after(tt,1)}{i}{5}{tsub}=[SelfInd{Ind_after(tt,2),Ind_after(tt,1)}{i}{5}{tsub}; reshape(temp_d(:,Ind_after(tt,2),Ind_after(tt,1),tsub)./NumSSDlet,1,NumElement)];
                        end
                    end
                end
                %% ====================================================================
                RM{index(j,2),index(j,1)}=L(index(j,2),index(j,1))*I_incident(j)*(I_after.*Wsub)'*M;% fluorescence emitted from current pixel
                TempSub(:,j,:)=L(index(j,2),index(j,1))*I_incident(j)*bsxfun(@times,I_after,M); % eqn(4) part 3
                xrfSub=xrfSub+RM{index(j,2),index(j,1)};
            end
            for j=1:size(index,1)
                temp=zeros(size(g(index(j,2),index(j,1),:)));
                %% ==================================================================== eqn(4) part 1
                if(j~=size(index,1))
                    for subj=j+1:size(index,1)
                        TempSub(:,j,:)=TempSub(:,j,:)- reshape(full(bsxfun(@times,reshape(MU_e(:,1,1).*L(index(j,2),index(j,1)),NumElement,1)...
                            ,RM{index(subj,2),index(subj,1)})),size(TempSub,1),1,size(TempSub,3));
                    end %%End loop for accumulating contribution from beam attenuation
                end
                %% ==================================================================== eqn(4) part 2
                if(~NoSelfAbsorption)
                    if(ismac)
                        TempIs=ismember(index(j,:),GlobalInd{n,i+1},'rows');
                    else
                        TempIs=ismember(index(j,:),GlobalInd{n,i+1},'rows','legacy');
                    end
                    if(~TempIs)
                        for i_sub=1:i
                            if(size(SelfInd{index(j,2),index(j,1)}{i_sub}{1},1)~=0)
                                TempLong=zeros(NumElement,numChannel);
                                for tsub=1:NumElement
                                    TempLong=TempLong+bsxfun(@times,(sum(bsxfun(@times,SelfInd{index(j,2),index(j,1)}{i_sub}{2}.*SelfInd{index(j,2),index(j,1)}{i_sub}{4}.*...
                                        SelfInd{index(j,2),index(j,1)}{i_sub}{3}(:,tsub),SelfInd{index(j,2),index(j,1)}{i_sub}{5}{tsub}),1))' ,M(tsub,:));
                                end
                                
                                TempSub(:,j,:)=TempSub(:,j,:)-reshape(TempLong,size(TempSub,1),1,size(TempSub,3));
                            end %% End loop for accumulating contribution from self absorption
                            clear TempLong
                        end
                    end
                end
                
                %% ====================================================================
                temp(1,1,:)=2*reshape(TempSub(:,j,:),NumElement,numChannel)*SigMa_XRF((nTau+1)*(n-1)+i)*(xrfSub-xrfData{n,i})';
                JacobSub1(index(j,2),index(j,1),:,:)=JacobSub1(index(j,2),index(j,1),:,:)+reshape(TempSub(:,j,:),1,1,NumElement,numChannel);
                g(index(j,2),index(j,1),:)=g(index(j,2),index(j,1),:)+alpha*temp;
                clear temp
            end
            if(LogScale)
                Rdis=e1'*(MU_XTM.*L)*e1; %% Discrete case
                sum_Tau=sum_Tau+alpha*SigMa_XRF((nTau+1)*(n-1)+i)*(xrfData{n,i}-xrfSub)*(xrfData{n,i}-xrfSub)'+beta*SigMa_XTM(i)*(Rdis-Mt(i))^2;
                g=g+2*beta*SigMa_XTM(i)*(Rdis-Mt(i)).*repmat(full(L),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]);
                JacobSub2=2*beta*SigMa_XTM(i)*repmat(full(L),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]);
            else
                Rdis=I0*exp(-e1'*(MU_XTM.*L)*e1);
                sum_Tau=sum_Tau+alpha*SigMa_XRF((nTau+1)*(n-1)+i)*(xrfData{n,i}-xrfSub)*(xrfData{n,i}-xrfSub)'+beta*SigMa_XTM(i)*(Rdis-Mt(i))^2;
                g=g-2*beta*SigMa_XTM(i)*Rdis*(Rdis-Mt(i)).*repmat(full(L),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]);
                JacobSub2=-2*beta*SigMa_XTM(i)*Rdis*repmat(full(L),[1,1,NumElement]).*repmat(MUe_XTM,[m(1),m(2),1]);
            end
            Jacob1=[Jacob1;reshape(JacobSub1,m(1)*m(2)*NumElement,numChannel)'];
            Jacob2=[Jacob2;JacobSub2(:)'];
        end
        
    end
    
    f=f+sum_Tau;
end
%%=============================================== Measure ill-conditioness
% figure,subplot(2,1,1),spy(Jacob1), subplot(2,1,2),spy(Jacob2);
% size(Jacob1)
% size(Jacob2)
% tol=eps;
% J1=rank(Jacob1);
% k1=cond(Jacob1);
% % [U1,S1,V1]=svd(Jacob1);
% J2=rank(Jacob2);
% k2=cond(Jacob2);
% % [U2,S2,V2]=svd(Jacob2);
% [U,S,V]=svd([Jacob1;Jacob2]);
% J3=rank([Jacob1;Jacob2]);
% fprintf('Cond1=%d, Rank1=%d, Cond2=%d, Rank2=%d, cond=%d, minS=%d, rank=%d \n',k1,J1,k2,J2,cond([Jacob1;Jacob2]),min(diag(S)),J3);
%%=====================================================================
g=g(:);
g=g(VarInd);
