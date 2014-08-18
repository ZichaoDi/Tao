function [f,g]=sfun_XRF_full2(W,xrfData,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,thetan,m,nTau)
global BeforeEmit  SSDlet dz omega NumSSDlet NoSelfAbsorption testind
f=0;
W=reshape(W,m(1),m(2),NumElement);
%%%%% =================== Attenuation Matrix at beam energy
% AbsorbScale=1;%1e-2;
% MU_e=MU_e.*AbsorbScale;
MU=zeros(m);
for i=1:m(1)
    for j=1:m(2)
        MU(i,j)=sum(reshape(W(i,j,:),NumElement,1).*reshape(MU_e(:,1,1),NumElement,1));
    end
end
%%%%% =================== Attenuation Matrix at flourescence energy (Corrected Attenuation)
MU_after=cell(NumElement,1);
for i=1:NumElement
    
    for t=1:m(1)
        for j=1:m(2)
            MU_after{i}(t,j)=sum(reshape(W(t,j,:),NumElement,1).*reshape(MU_e(:,1,i+1),NumElement,1));
        end
    end
end
%%%%% ====================================================================
g=zeros(m(1),m(2),NumElement);
SelfInd=cell(m(1),m(2));
for im=1:m(1)
    for jm=1:m(2)
        SelfInd{im,jm}=cell(nTau+1,1);
        for i=1:nTau+1
            SelfInd{im,jm}{i}=cell(5,1);
        end
    end
end
clear i j t
for n=1:length(thetan)
    theta=thetan(n)/180*pi;
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    SSDknot=SSDlet*TransMatrix;
    sum_Tau=0;
    for i=1:nTau+1
        index=GlobalInd{n,i};
        if(~isempty(index))
            L=Ltol{n,i};
            RM=cell(m(1),m(2));
            xrfSub=zeros(1,numChannel);
            TempSub=zeros(NumElement,size(index,1),numChannel);
            I_incident=[];
            temp_d=zeros(NumElement,m(1),m(2));
            for j=1:size(index,1)
                CurrentCellCenter=[(index(j,1)-1/2)*dz-abs(omega(1)),(index(j,2)-1/2)*dz-abs(omega(3))];
                
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
                    for tsub=1:NumElement
                        for SSDi=1:NumSSDlet
                            temp_after=0;
                            beta=angle(SSDknot(SSDi,1)-CurrentCellCenter(1)+(SSDknot(SSDi,2)-CurrentCellCenter(2))*sqrt(-1));
                            if(beta>=0 & beta<=pi/2)
                                xbox=[CurrentCellCenter(1) CurrentCellCenter(1) omega(2) omega(2) CurrentCellCenter(1)];
                                ybox=[CurrentCellCenter(2) omega(4) omega(4) CurrentCellCenter(2) CurrentCellCenter(2)];
                            elseif(beta >pi/2 & beta<=pi)
                                xbox=[omega(1) omega(1) CurrentCellCenter(1) CurrentCellCenter(1) omega(1)];
                                ybox=[ CurrentCellCenter(2) omega(4) omega(4) CurrentCellCenter(2) CurrentCellCenter(2)];
                            elseif(beta >=-pi/2 & beta<0)
                                xbox=[CurrentCellCenter(1) CurrentCellCenter(1) omega(2) omega(2) CurrentCellCenter(1)];
                                ybox=[omega(3) CurrentCellCenter(2) CurrentCellCenter(2) omega(3) omega(3)];
                            elseif(beta>=-pi & beta<-pi/2)
                                xbox=[omega(1) omega(1) CurrentCellCenter(1) CurrentCellCenter(1) omega(1)];
                                ybox=[omega(3) CurrentCellCenter(2) CurrentCellCenter(2) omega(3) omega(3)];
                            end
                            if(beta<0)
                                beta=beta+2*pi;
                            end
                            [index_after,Lvec_after]=IntersectionSet(CurrentCellCenter,SSDknot(SSDi,:),xbox,ybox,beta);
                            [index_after,otherInd]=setdiff(index_after,index(j,:),'rows');
                            Ind_after=[Ind_after;index_after];
                            Lvec_after=Lvec_after(otherInd);
                            for t_after=1:size(index_after,1)
                                temp_after=temp_after+Lvec_after(t_after)*MU_after{tsub}(index_after(t_after,2),index_after(t_after,1));
                            end %% Attenuation of Flourescent energy emitted from current pixel
                            for si=1:size(index_after,1)
                                temp_d(:,index_after(si,2),index_after(si,1))=temp_d(:,index_after(si,2),index_after(si,1))+exp(-temp_after)*reshape(MU_e(:,1,tsub+1),NumElement,1).*Lvec_after(si);
                                if(~ismember(index(j,end:-1:1),SelfInd{index_after(si,2),index_after(si,1)}{i}{1},'rows'))
                                    SelfInd{index_after(si,2),index_after(si,1)}{i}{1}=[SelfInd{index_after(si,2),index_after(si,1)}{i}{1};index(j,end:-1:1)];%% assign downstream index to pixel(oppsite)
                                    SelfInd{index_after(si,2),index_after(si,1)}{i}{2}=[SelfInd{index_after(si,2),index_after(si,1)}{i}{2};L(index(j,2),index(j,1))];
                                    SelfInd{index_after(si,2),index_after(si,1)}{i}{3}=[SelfInd{index_after(si,2),index_after(si,1)}{i}{3};reshape(W(index(j,2),index(j,1),:),1,NumElement)];
                                    SelfInd{index_after(si,2),index_after(si,1)}{i}{4}=[SelfInd{index_after(si,2),index_after(si,1)}{i}{4};I_incident(j)];
                                end
                            end
                            I_after(tsub)=I_after(tsub)+exp(-temp_after);
                        end %% End loop for each SSD detector let
                        Ind_after=unique(Ind_after,'rows');
                        I_after(tsub)=I_after(tsub)/NumSSDlet;
                    end %% End loop for existing fluorescence energy from current pixel
                    for tt=1:size(Ind_after,1)
                        SelfInd{Ind_after(tt,2),Ind_after(tt,1)}{i}{5}=[SelfInd{Ind_after(tt,2),Ind_after(tt,1)}{i}{5};reshape(temp_d(:,Ind_after(tt,2),Ind_after(tt,1))./NumSSDlet,1,NumElement)];
                    end
                end
                
                %% ====================================================================
                RM{index(j,2),index(j,1)}=L(index(j,2),index(j,1))*I_incident(j)*(I_after.*Wsub)'*M;% fluorescence emitted from current pixel
                TempSub(:,j,:)=L(index(j,2),index(j,1))*I_incident(j)*bsxfun(@times,I_after,M);%repmat(I_after,1,size(M,2))).*M; % eqn(4) part 3
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
                    if(~ismember(index(j,:),GlobalInd{n,i+1},'rows'))
                        for i_sub=1:i
                            for diffj_i=1:size(SelfInd{index(j,2),index(j,1)}{i_sub}{1},1)
                                TempLong=sum(SelfInd{index(j,2),index(j,1)}{i_sub}{2}(diffj_i)*SelfInd{index(j,2),index(j,1)}{i_sub}{4}(diffj_i).*...
                                    bsxfun(@times,SelfInd{index(j,2),index(j,1)}{i_sub}{3}(diffj_i,:)',M),1);
                                TempSub(:,j,:)=TempSub(:,j,:)-reshape(full(bsxfun(@times,SelfInd{index(j,2),index(j,1)}{i_sub}{5}(diffj_i,:)',TempLong)),size(TempSub,1),1,size(TempSub,3));
                                
                                clear TempLong
                            end %% End loop for accumulating contribution from self absorption
                        end
                    end
                end
                %% ====================================================================
                temp(1,1,:)=2*reshape(TempSub(:,j,:),NumElement,numChannel)*(xrfSub-xrfData{n,i})';
                g(index(j,2),index(j,1),:)=g(index(j,2),index(j,1),:)+temp;
                clear temp
            end
            sum_Tau=sum_Tau+(xrfData{n,i}-xrfSub)*(xrfData{n,i}-xrfSub)';
        end
    end
    
    f=f+sum_Tau;
end

g=g(:);
