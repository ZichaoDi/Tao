%%%Simulate XRF of a given object with predifined detector and beam
% function XRF=SimulateXRF(W,MU,BindingEenergy,M,thetan,DetChannel, numChannel, nTau, DetKnot0, SourceKnot0);
global plotSpecSingle BeforeEmit plotTravel SSDlet NoSelfAbsorption
global fig2  fig5 finalfig
% close all;
plotTravel=0; % If plot the intersection of beam with object
plotSpec = 0; % Do you want to see the spectra? If so plotSpec = 1
plotUnit=0;
plotSpecSingle=0;
NoSelfAbsorption=1;
PlotObject=1;
% global m x y omega NumElement plotTravel plotSpec plotWhole dz
startup;
more off;

Define_Detector_Beam_Gaussian; %% provide the beam source and Detectorlet
UnitSpectrumSherman_Gaussian; %% Produce BindingEnergy M
DefineObject_Gaussian; %% Produce W, MU
thetan=[1:40:180];%linspace(1,180,25);%[1:60:180];% Projection Angles
%%%%%%%==============================================================
if plotTravel
    fig2=[];  fig5=[];
end
%%%=========Start Simulation
Energy=BindingEnergy;
e1=ones(m(1),1);
XRF=cell(length(thetan),nTau+1);
DisR=zeros(nTau+1,length(thetan));
Ltol=cell(length(thetan),nTau+1);
GlobalInd=cell(length(thetan),nTau+1);
LocalInd=cell(length(thetan),nTau+1,m(1),m(2),NumSSDlet);
L_after=cell(length(thetan),nTau+1,m(1),m(2),NumSSDlet);
RMlocal=zeros(m(1),m(2),numChannel); %% assign all the contributions from seperate beam to each pixel
fprintf(1,'====== Transmission Detector Resolution is %d\n',nTau);
fprintf(1,'====== Fluorescence Detector Resolution is %d\n',numChannel);
for n=1:length(thetan)
    
    theta=thetan(n)/180*pi;
    fprintf(1,'====== Angle Number  %d of %d: %d\n',n,length(thetan),thetan(n));
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    DetKnot=DetKnot0*TransMatrix;
    SourceKnot=SourceKnot0*TransMatrix;
    SSDknot=SSDlet*TransMatrix;
    %%%%%============================================================
%         fig=figure('name','grids & beam line');
%         plotGrid(xc,omega,m); hold on;
%         for i=1:nTau+1
%             fig3=plot([SourceKnot(i,1),DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)],'r.-');hold on;
%         end
%         pause;
    %%%%%%%===============================================
    Rdis=zeros(nTau+1,1);
    for i=1:nTau+1 %%%%%%%%%================================================================
        % Initialize
        XRF{n,i} = zeros(numChannel,1);
        
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        BeforeEmit=1;
        [index,Lvec]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta);
        L=sparse(zeros(m(1),m(2)));
        %============================= Plot Grid and Current Light Beam
        
        if(plotSpec)
            finalfig=figure('name',sprintf('XRF at beam %d with angle %d',i,thetan(n)));
            subplot(1,2,1);
            plotGrid(xc,omega,m);
            hold on;
            plot(DetKnot(:,1),DetKnot(:,2),'k+-',SourceKnot(:,1),SourceKnot(:,2),'m+-',SSDknot(:,1),SSDknot(:,2),'g+-','LineWidth',0.5)
            axis equal
            set(gcf,'Units','normalized')
            set(gca,'Units','normalized')
            ax = axis;
            ap = get(gca,'Position');
            xp = ([SourceKnot(i,1),DetKnot(i,1)]-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
            yp = ([SourceKnot(i,2),DetKnot(i,2)]-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
            ah=annotation('arrow',xp,yp,'Color','r','LineStyle','--');
        end
        %=================================================================
        %%%%%%%%================================================================
        GlobalInd{n,i}=index;
        RM=cell(m(1),m(2));
        xrfSub=zeros(1,numChannel);
        for j=1:size(index,1)
            CurrentCellCenter=[(index(j,1)-1/2)*dz-abs(omega(1)),(index(j,2)-1/2)*dz-abs(omega(3))];
            L(index(j,2),index(j,1))=Lvec(j);
            if(j==1)
                I_incident=1;
                temp_sum=0;
            else
                temp_sum=temp_sum+Lvec(j-1)*MU(index(j-1,2),index(j-1,1));
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
                    [index_after,Lvec_after]=IntersectionSet(CurrentCellCenter,SSDknot(SSDi,:),xbox,ybox,beta);
                    [index_after,otherInd]=setdiff(index_after,index(j,:),'rows');
                    Lvec_after=Lvec_after(otherInd);
                    LocalInd{n,i,index(j,2),index(j,1),SSDi}=index_after;
                    L_after{n,i,index(j,2),index(j,1),SSDi}=Lvec_after;
                    for tsub=1:NumElement
                        for t_after=1:size(index_after,1)
                            temp_after=temp_after+Lvec_after(t_after)*MU_after{tsub}(index_after(t_after,2),index_after(t_after,1));
                        end %% Attenuation of Flourescent energy emitted from current pixel
                        I_after(tsub)=I_after(tsub)+exp(-temp_after)/NumSSDlet;
                    end %% End loop for each SSD detector let
                end %% End loop for existing fluorescence energy from current pixel
            end
            %% ====================================================================
            RM{index(j,2),index(j,1)}=Lvec(j)*I_incident*(I_after.*Wsub)'*M;% fluorescence emitted from current pixel
            temp=zeros(size(RMlocal(index(j,2),index(j,1),:)));
            temp(1,1,:)=Lvec(j)*I_incident*(I_after.*Wsub)'*M;
            RMlocal(index(j,2),index(j,1),:)=RMlocal(index(j,2),index(j,1),:)+temp;
            xrfSub=xrfSub+RM{index(j,2),index(j,1)};
        end
        Ltol{n,i}=L;
        Rdis(i)=I0*exp(-e1'*(MU.*L)*e1); %% Discrete case
        XRF{n,i}=xrfSub;
        if(plotSpec)
            figure(finalfig)
            subplot(1,2,2);
            plot(DetChannel,XRF{n,i},'r-')
            xlabel('Energy Channel','fontsize',12); ylabel('Intensity','fontsize',12)
            pause;
        end
        
    end
    DisR(:,n)=Rdis';
end

