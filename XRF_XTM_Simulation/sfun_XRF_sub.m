function [f,g]=sfun_XRF_sub(W,xrfData,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,thetan,m,nTau)
global BeforeEmit  SSDlet dz omega NumSSDlet NoSelfAbsorption AbsorbScale
global testind
f=0;
W=reshape(W,m(1),m(2),NumElement);
%%%%% =================== Attenuation Matrix at beam energy
% AbsorbScale=1e0;
MU_e=MU_e.*AbsorbScale;
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
g=zeros(m(1),m(2),NumElement,numChannel);
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
            temp_d=ones(NumElement,size(index,1));
            SelfInd=cell(size(index,1),1);
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
                        Lvec_after=Lvec_after(otherInd);
                        [~,tempind1,tempind2]=intersect(index,index_after,'rows');
                        for t_after=1:size(index_after,1)
                            temp_after=temp_after+Lvec_after(t_after)*MU_after{tsub}(index_after(t_after,2),index_after(t_after,1));
                        end %% Attenuation of Flourescent energy emitted from current pixel
                        for si=1:length(tempind1)
                        temp_d(:,tempind1(si))=temp_d(:,tempind1(si))+exp(-temp_after)*reshape(MU_e(:,1,tsub+1),NumElement,1)...
                            .*Lvec_after(tempind2(si))/NumSSDlet;
                        SelfInd{tempind1(si)}=[SelfInd{tempind1(si)};j];%% assign downstream index to pixel(oppsite)
                        
                        end
                        I_after(tsub)=I_after(tsub)+exp(-temp_after);
                    end %% End loop for each SSD detector let
                    I_after(tsub)=I_after(tsub)/NumSSDlet;
                end %% End loop for existing fluorescence energy from current pixel
                end

                %% ====================================================================
                RM{index(j,2),index(j,1)}=L(index(j,2),index(j,1))*I_incident(j)*(I_after.*Wsub)'*M;% fluorescence emitted from current pixel
                TempSub(:,j,:)=(L(index(j,2),index(j,1))*I_incident(j)*repmat(I_after,1,size(M,2))).*M; % eqn(4) part 3
                xrfSub=xrfSub+RM{index(j,2),index(j,1)};
            end
            for j=1:size(index,1)
                temp=zeros(1,1,NumElement,numChannel);
                %% ==================================================================== eqn(4) part 1
                if(j~=size(index,1))
                    for subj=j+1:size(index,1)
                        TempSub(:,j,:)=TempSub(:,j,:)- reshape(full(repmat(reshape(MU_e(:,1,1).*L(index(j,2),index(j,1)),NumElement,1),1,size(M,2))...
                            .*repmat(RM{index(subj,2),index(subj,1)},NumElement,1)),size(TempSub,1),1,size(TempSub,3));
                    end %%End loop for accumulating contribution from beam attenuation
                end
                
                %% ==================================================================== eqn(4) part 2
                if(~NoSelfAbsorption)
                    diffj=unique(SelfInd{j});
                    for diffj_i=1:length(diffj)                            
                        
                        TempLong=sum(L(index(diffj(diffj_i),2),index(diffj(diffj_i),1))*I_incident(diffj(diffj_i)).*...
                            repmat(reshape(W(index(diffj(diffj_i),2),index(diffj(diffj_i),1),:),NumElement,1),1,size(M,2)).*M,1);
                        TempSub(:,j,:)=TempSub(:,j,:)-reshape(full(repmat(temp_d(:,diffj(diffj_i)),1,size(M,2)).*repmat(TempLong,NumElement,1)),size(TempSub,1),1,size(TempSub,3));
                     clear TempLong   
                    end %% End loop for accumulating contribution from self absorption
                end
                %% ====================================================================
                
                temp(1,1,:,:)=reshape(TempSub(:,j,:),NumElement,numChannel);
                g(index(j,2),index(j,1),:,:)=g(index(j,2),index(j,1),:,:)+temp;
                clear temp
            end
            sum_Tau=sum_Tau+xrfSub;
        end
    end
   
    f=f+sum_Tau;
end
g=reshape(g,m(1)*m(2)*NumElement,numChannel);
