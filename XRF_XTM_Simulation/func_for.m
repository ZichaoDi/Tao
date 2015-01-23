function f=func_for(W,xrfData,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,LocalInd,L_after,thetan,m,nTau)
global BeforeEmit  NumSSDlet NoSelfAbsorption
global SigMa_XRF
f=0;
W=reshape(W,m(1),m(2),NumElement);
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,1,NumElement);
MU=sum(W.*repmat(MUe,[m(1),m(2),1]),3);
%%%%% =================== Attenuation Matrix at flourescence energy (Corrected Attenuation)
MU_after=cell(NumElement,1);
for i=1:NumElement
    MU_after{i}=sum(W.*repmat(reshape(MU_e(:,1,i+1),1,1,NumElement),[m(1),m(2),1]),3);
end
%%%%% ====================================================================
for n=1:length(thetan)
    sum_Tau=0;
    for i=1:nTau+1
        index=GlobalInd{n,i};
        if(~isempty(index))
            L=Ltol{n,i};
            RM=cell(m(1),m(2));
            xrfSub=zeros(1,numChannel);
            I_incident=zeros(size(index,1),1);
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
                            if(~isempty(Lvec_after))
                            temp_after=sum(Lvec_after.*reshape(MU_after{tsub}(LinearInd),size(Lvec_after))); %% Attenuation of Flourescent energy emitted from current pixel
                            else
                                temp_after=0;
                            end
                            I_after(tsub)=I_after(tsub)+exp(-temp_after)/NumSSDlet;
                        end %% End loop for existing fluorescence energy from current pixel
                    end %% End loop for each SSD detector let
                end
                %% ====================================================================
                RM{index(j,2),index(j,1)}=L(index(j,2),index(j,1))*I_incident(j)*(I_after.*Wsub)'*M;% fluorescence emitted from current pixel
                xrfSub=xrfSub+RM{index(j,2),index(j,1)};
            end
        sum_Tau=sum_Tau+SigMa_XRF((nTau+1)*(n-1)+i)*(xrfData{n,i}-xrfSub)*(xrfData{n,i}-xrfSub)';
        end
    end
end

f=f+sum_Tau;
end

