%%%Simulate XRF of a given object with predifined detector and beam
close all;
load Phantom3_Young
%Phantom3_Young=Phantom3_Young(100:120,100:120,100:120,:);
onlyXRF=1;
more off;
current_n=size(Phantom3_Young,1);
xrf_roi=0;
NumElement=4;
numThetan=10;
DefineGeometry;
DefineObject_Gaussian; % Produce W, MU_XTM
%%-------------------------------------------------------
dataset={'channel_names', 'channel_units','scaler_names', 'scaler_units', 'scalers','XRF_roi','mca_arr'};
load formatH5_roi
%%---------------------------------------------------------------------------
TotSlice=current_n;
XRF=zeros(numThetan,nTau+1,numChannel);
%X1=cell(numThetan,1);
%X2=cell(numThetan,1);
disp('Setup Geometry')
for n=1:length(thetan)
    n
    XRF3=zeros(nTau+1,TotSlice,numChannel);
    XRF_roi=zeros(nTau+1,TotSlice,NumElement);
    for slice=1:TotSlice
        W=zeros(m(1),m(2),NumElement);
        for NE=1:NumElement
            data=Phantom3_Young(:,:,:,NE+1);
            W(:,:,NE)=data(:,:,slice);
        end
        
        %%%%% =================== Attenuation Matrix at beam energy
        MU_XTM=sum(W.*repmat(reshape(MU_e(:,1,1),1,1,NumElement),[m(1),m(2),1]),3);
        %%%%% =================== Attenuation Matrix at flourescence energy (Corrected Attenuation)
        MU_after=cell(NumElement,1);
        for i=1:NumElement
            MU_after{i}=sum(W.*repmat(reshape(MU_e(:,1,i+1),1,1,NumElement),[m(1),m(2),1]),3);
        end
        %%%%% ====================================================================
        %%---------------------------------------------------------------------------
        for i=1:nTau+1 %%%%%%%%%================================================================
            BeforeEmit=1;
            %=================================================================
            index=ID{n,i};
            Lvec=LD{n,i};
            xrfSub=zeros(1,numChannel);
            %%%%%%%%================================================================
            for j=1:size(index,1)
                if(j==1)
                    I_incident=1;
                    temp_sum=0;
                elseif(j>1 & j<size(index,1))
                    temp_sum=temp_sum+Lvec(j-1)*MU_XTM(index(j-1,2),index(j-1,1));
                    I_incident=exp(-temp_sum);
                else
                    temp_sum=temp_sum+Lvec(j-1)*MU_XTM(index(j-1,2),index(j-1,1));
                    I_incident=exp(-temp_sum);
                end
                %% ===========================================================
                Wsub=reshape(W(index(j,2),index(j,1),:),[NumElement,1]);
                %% Self-absorption
                I_after=ones(NumElement,1);
                for SSDi=1:NumSSDlet
                    temp_after=0;
                    Lvec_after=LA{n,i,index(j,2),index(j,1),SSDi};
                    LinearInd=LI{n,i,index(j,2),index(j,1),SSDi};
                    for tsub=1:NumElement
                        if(~isempty(Lvec_after))
                            temp_after=sum(Lvec_after.*reshape(MU_after{tsub}(LinearInd),size(Lvec_after))); %% Attenuation of Flourescent energy emitted from current pixel
                        else
                            temp_after=0;
                        end
                        I_after(tsub)=I_after(tsub)+exp(-temp_after)/NumSSDlet;
                    end %% End loop for each SSD detector let
                end %% End loop for existing fluorescence energy from current pixel
                %% ====================================================================
                xrfSub=xrfSub+Lvec(j)*I_incident*(I_after.*Wsub)'*M;% fluorescence emitted from current pixel
            end
            XRF(n,i,:)=xrfSub;
            for whichElement=1:NumElement
                XRF_roi(i,slice,whichElement)=sum(XRF(n,i,EmitEner(whichElement,:)),3);
            end
        end
        XRF3(:,slice,:)=reshape(XRF(n,:,:),[nTau+1,1,numChannel]);
    end
    formatH5_roi{length(dataset)}=XRF3;
    formatH5_roi{length(dataset)-1}=XRF_roi;
    for ii=1:length(dataset)-1
        if(ii==1)
            hdf5write(['YoungPhantom_256_4_',num2str(n),'.h5'],['/MAPS/',num2str(dataset{ii})],formatH5_roi{ii});
        else
           hdf5write(['YoungPhantom_256_4_',num2str(n),'.h5'],['/MAPS/',num2str(dataset{ii})],formatH5_roi{ii},'WriteMode', 'append');
        end
    end
    hdf5write(['YoungPhantom_256_4_',num2str(n),'.h5'],'/MAPS/mca_arr',formatH5_roi{end},'WriteMode', 'append');
   % X1{n}=XRF3;
%X2{n}=XRF_roi;

%save -v7.3 X1 X1
%save -v7.3 X2 X2
end





