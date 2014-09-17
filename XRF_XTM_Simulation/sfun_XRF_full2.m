function [f,g]=sfun_XRF_full2(W,xrfData,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,LocalInd,L_after,thetan,m,nTau)
global BeforeEmit  NumSSDlet NoSelfAbsorption testind
global SigMa_XRF
f=0;
W=reshape(W,m(1),m(2),NumElement);
%%%%% =================== Attenuation Matrix at beam energy
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,1,NumElement);
MU=sum(W.*repmat(MUe,[m(1),m(2),1]),3);
%%%%% =================== Attenuation Matrix at flourescence energy (Corrected Attenuation)
MU_after=cell(NumElement,1);
for i=1:NumElement
MU_after{i}=sum(W.*repmat(reshape(MU_e(:,1,i+1),1,1,NumElement),[m(1),m(2),1]),3);
end
%%%%% ====================================================================

g=zeros(m(1),m(2),NumElement);
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
for n=1:length(thetan)
    SelfInd=SelfInd1;
    sum_Tau=0;
    for i=1:nTau+1
        index=GlobalInd{n,i};
        if(~isempty(index))
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
                if(NoSelfAbsorption)
                    NumSSDlet=1;%% Turn off self-absorption
                else
                    I_after=0*I_after;
                    
                    for SSDi=1:NumSSDlet
                        index_after=LocalInd{n,i,index(j,2),index(j,1),SSDi};
                        Lvec_after=L_after{n,i,index(j,2),index(j,1),SSDi};
                        LinearInd=sub2ind([m(1),m(2)],index_after(:,2),index_after(:,1));
                        for tsub=1:NumElement
                            temp_after=sum(Lvec_after.*MU_after{tsub}(LinearInd)); %% Attenuation of Flourescent energy emitted from current pixel
%                             LinearInd1=sub2ind([NumElement,m(1),m(2),NumElement],repmat(1:NumElement,1,length(index_after(:,2))),repmat(index_after(:,2)',1,NumElement),...
%                                 repmat(index_after(:,1)',1,NumElement),tsub*ones(1,NumElement*length(index_after(:,2))));
%                             temp_d_sub=bsxfun(@times,reshape(MU_e(:,1,tsub+1),NumElement,1),Lvec_after');
%                             temp_d(LinearInd1)=exp(-temp_after)*temp_d_sub(:);
                            for si=1:size(index_after,1)
                                temp_d(:,index_after(si,2),index_after(si,1),tsub)=temp_d(:,index_after(si,2),index_after(si,1),tsub)+exp(-temp_after)*reshape(MU_e(:,1,tsub+1),NumElement,1).*Lvec_after(si);
                                if(~ismember(index(j,end:-1:1),SelfInd{index_after(si,2),index_after(si,1)}{i}{1},'rows'))
                                    SelfInd{index_after(si,2),index_after(si,1)}{i}{1}=[SelfInd{index_after(si,2),index_after(si,1)}{i}{1};index(j,end:-1:1)];%% assign downstream index to pixel(oppsite)
                                    SelfInd{index_after(si,2),index_after(si,1)}{i}{2}=[SelfInd{index_after(si,2),index_after(si,1)}{i}{2};L(index(j,2),index(j,1))];
                                    SelfInd{index_after(si,2),index_after(si,1)}{i}{3}=[SelfInd{index_after(si,2),index_after(si,1)}{i}{3};reshape(W(index(j,2),index(j,1),:),1,NumElement)];
                                    SelfInd{index_after(si,2),index_after(si,1)}{i}{4}=[SelfInd{index_after(si,2),index_after(si,1)}{i}{4};I_incident(j)];
                                end
                            end
                            I_after(tsub)=I_after(tsub)+exp(-temp_after)/NumSSDlet;
                        end %% End loop for existing fluorescence energy from current pixel
                    end %% End loop for each SSD detector let
                    Ind_after=unique(cat(1,LocalInd{n,i,index(j,2),index(j,1),:}),'rows');
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
                TempSub(:,j,:)=L(index(j,2),index(j,1))*I_incident(j)*bsxfun(@times,I_after,M); % eqn(5) part 3
                xrfSub=xrfSub+RM{index(j,2),index(j,1)};
            end
            for j=1:size(index,1)
                temp=zeros(size(g(index(j,2),index(j,1),:)));
                %% ==================================================================== eqn(5) part 1
                if(j~=size(index,1))
                    for subj=j+1:size(index,1)
                        TempSub(:,j,:)=TempSub(:,j,:)- reshape(full(bsxfun(@times,reshape(MU_e(:,1,1).*L(index(j,2),index(j,1)),NumElement,1)...
                            ,RM{index(subj,2),index(subj,1)})),size(TempSub,1),1,size(TempSub,3));
                    end %%End loop for accumulating contribution from beam attenuation
                end
                %% ==================================================================== eqn(5) part 2
                if(~NoSelfAbsorption)
                    if(~ismember(index(j,:),GlobalInd{n,i+1},'rows'))
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
                temp(1,1,:)=2*reshape(TempSub(:,j,:),NumElement,numChannel)*SigMa_XRF((nTau+1)*(n-1)+i)*(xrfSub-xrfData{n,i})';%
                g(index(j,2),index(j,1),:)=g(index(j,2),index(j,1),:)+temp;
                clear temp
            end
            sum_Tau=sum_Tau+SigMa_XRF((nTau+1)*(n-1)+i)*(xrfData{n,i}-xrfSub)*(xrfData{n,i}-xrfSub)';%
        end
    end
    
    f=f+sum_Tau;
end

g=g(:);
