%%%Simulate XRF of a given object with predifined detector and beam
% function XRF=SimulateXRF(W,MU,BindingEenergy,M,thetan,DetChannel, numChannel, nTau, DetKnot0, SourceKnot0);
global plotSpecSingle BeforeEmit plotTravel
global fig2  fig5 finalfig slice onlyXRF

Tomo_startup;
close all;
load Phantom3_Young
onlyXRF=1;
PlotObject=0;
plotSpec = 0; % Do you want to see the spectra? If so plotSpec = 1
plotTravel=0; % If plot the intersection of beam with object
plotUnit=0;
plotElement=0;
plotSpecSingle=0;
NoSelfAbsorption=0;
more off;
current_n=1;
numThetan=1;
XRF3=zeros(nTau+1,current_n,numChannel);
xrf_roi=0;
Define_Detector_Beam_Gaussian; %% provide the beam source and Detectorlet
DefineObject_Gaussian; % Produce W, MU_XTM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XRF=cell(length(thetan),nTau+1);
for n=1:length(thetan)
    theta=thetan(n)/180*pi;
    fprintf(1,'====== Angle Number  %d of %d: %d\n',n,length(thetan),thetan(n));
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    DetKnot=DetKnot0*TransMatrix;
    SourceKnot=SourceKnot0*TransMatrix;
    SSDknot=SSDlet*TransMatrix;
    for i=1:nTau+1 %%%%%%%%%================================================================
        XRF{n,i} = zeros(numChannel,1);   
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        BeforeEmit=1;
        %=================================================================
        [index,Lvec,linearInd]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta);
        xrfSub=zeros(1,numChannel);
        %%%%%%%%================================================================
        for j=1:size(index,1)
            CurrentCellCenter=[(index(j,1)-1/2)*dz(1)-abs(omega(1)),(index(j,2)-1/2)*dz(2)-abs(omega(3))];
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
            
            BeforeEmit=0;
            I_after=ones(NumElement,1);
            if(NoSelfAbsorption)
                NumSSDlet=1;%% Turn off self-absorption
            else
                I_after=0*I_after;
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
                    [index_after,Lvec_after,linearInd_after]=IntersectionSet(CurrentCellCenter,SSDknot(SSDi,:),xbox,ybox,beta);
                    [index_after,otherInd]=setdiff(index_after,index(j,:),'rows');
                    Lvec_after=Lvec_after(otherInd');
                    LinearInd=sub2ind(m,index_after(:,2),index_after(:,1));
                    
                    
  for slice=1:current_n
MU_after=cell(NumElement,1);
for i=1:NumElement
    MU_after{i}=sum(W.*repmat(reshape(MU_e(:,1,i+1),1,1,NumElement),[m(1),m(2),1]),3);
end
W=zeros(current_n,current_n,NumElement);
        for i=1:NumElement
            data=Phantom3_Young(:,:,:,i+1);
            W(:,:,i)=data(1:current_n,1:current_n,slice);
        end                  
                    
                    for tsub=1:NumElement
                        if(~isempty(Lvec_after))
                            temp_after=sum(Lvec_after.*reshape(MU_after{tsub}(LinearInd),size(Lvec_after))); %% Attenuation of Flourescent energy emitted from current pixel
                        else
                            temp_after=0;
                        end
                        I_after(tsub)=I_after(tsub)+exp(-temp_after)/NumSSDlet;
                    end %% End loop for each SSD detector let
                end %% End loop for existing fluorescence energy from current pixel
            end
            %% ====================================================================
            xrfSub=xrfSub+Lvec(j)*I_incident*(I_after.*Wsub)'*M;% fluorescence emitted from current pixel
            
        end
        XRF{n,i}=xrfSub;
    end
end

XRF3(1:nTau+1,slice,1:numChannel)=reshape(cell2mat(XRF'),[nTau+1,1,numChannel]);
end


dataset={'channel_names', 'channel_units','scaler_names', 'scaler_units', 'scalers','XRF_roi','mca_arr'};
% formatH5=cell(length(dataset),1);
% for i=1:length(dataset)-1
% formatH5{i}=h5read('/Users/Wendydi/Documents/MATLAB/APSdata/2xfm1211_14/2xfm_0307.h5',['/MAPS/',num2str(dataset{i})]);
% end
load formatH5_roi
formatH5_roi{length(dataset)}=XRF3;
% for i=1:length(dataset)
% if(i==1)
% hdf5write('myXRF.h5',['/MAPS/',num2str(dataset{i})],formatH5{i});
% else
% hdf5write('myXRF.h5','/MAPS/mca_arr',formatH5_roi{end},'WriteMode', 'append');
hdf5write('myXRF1.h5','/MAPS/mca_arr',formatH5_roi{end});
% end
% end

