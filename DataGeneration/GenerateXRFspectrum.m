%%%Simulate XRF of a given object with predifined detector and beam
close all;
load Phantom3_Young
onlyXRF=1;
more off;
current_n=size(Phantom3_Young,1);
xrf_roi=0;
NumElement=4;
Define_Detector_Beam_Gaussian; %% provide the beam source and Detectorlet
thetan=linspace(0,180,21);%mod(thetan+360,360);% Projection Angles, has to be positive.
%thetan=thetan(15:end);
numThetan=length(thetan);
DefineObject_Gaussian; % Produce W, MU_XTM
for n=1:numThetan
    DefineGeometry;
    %%-------------------------------------------------------
    dataset={'channel_names', 'channel_units','scaler_names', 'scaler_units', 'scalers','XRF_roi','mca_arr'};
    load formatH5_roi
    %%---------------------------------------------------------------------------
    TotSlice=current_n;
    disp('Setup Geometry')
    XRF=zeros(nTau+1,numChannel);
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
            %=================================================================
            index=ID{i};
            Lvec=LD{i};
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
                %% Modeling Self-absorption
                I_after=ones(NumElement,1);
                for SSDi=1:NumSSDlet
                    temp_after=0;
                    Lvec_after=LA{i,index(j,2),index(j,1),SSDi};
                    LinearInd=LI{i,index(j,2),index(j,1),SSDi};
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
            XRF(i,:)=xrfSub;
            for whichElement=1:NumElement
                XRF_roi(i,slice,whichElement)=sum(XRF(i,EmitEner(whichElement,:)),2);
            end
        end
        XRF3(:,slice,:)=reshape(XRF(:,:),[nTau+1,1,numChannel]);
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
end




