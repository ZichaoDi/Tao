%%%Simulate XRF of a given object with predifined detector and beam
% function XRF=SimulateXRF(W,MU,BindingEenergy,M,thetan,DetChannel, numChannel, nTau, DetKnot0, SourceKnot0);
global plotSpecSingle BeforeEmit plotTravel SSDknot fig
global fig1 fig2 fig4 fig5 fig6
close all;
plotTravel=1; % If plot the intersection of beam with object
plotSpec = 1; % Do you want to see the spectra? If so plotSpec = 1
plotWhole =0;
plotUnit=0;
plotSpecSingle=0;
% global m x y omega NumElement plotTravel plotSpec plotWhole dz
startup;
more off;
DefineObject; %% Produce W, MU
Define_Detector_Beam; %% provide the beam source and Detectorlet
% UnitSpectrum; %% Produce BindingEnergy M
% UnitSpectrumSherman
thetan=[30 120 210 300];% Projection Angles
thetan=300;
%%%%%%%==============================================================
if plotTravel
    figure(fig);
    fig1=[]; fig2=[]; fig3=[]; fig4=[]; fig5=[]; fig6=[];
end
%%%=========Start Simulation
Energy=BindingEnergy;

XRF=cell(length(thetan),nTau+1);

Ltol=cell(length(thetan),nTau+1);
GlobalInd=cell(length(thetan),nTau+1);
RMlocal=zeros(m(1),m(2),NumElement);
fprintf(1,'====== Detector Resolution is %d\n',nTau);
for n=1:length(thetan)
    fprintf(1,'====== Angle Number  %d of %d\n',n,length(thetan));
    
    theta=thetan(n)/180*pi;
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    DetKnot=DetKnot0*TransMatrix;
    SourceKnot=SourceKnot0*TransMatrix;
    SSDknot=SSDlet*TransMatrix;
    if plotWhole
        figure(fig);
        set(fig3,'visible','off');
        for i=1:nTau+1
            fig3=plot([SourceKnot(i,1),DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)],'r.-',SourceKnot(i,1),SourceKnot(i,2),'bo');hold on;
        end
    end
    
    
    for i=1:nTau+1 %%%%%%%%%================================================================
        % Initialize
        XRF{n,i} = zeros(numChannel,1);
        
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        BeforeEmit=1;
        [index,Lvec]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta);
        if(isempty(index))
            L=[];
        else
            %%%%%%%%%================================================================
            GlobalInd{n,i}=index;            
            
            L=sparse(zeros(m(1),m(2)));
            RM=cell(m(1),m(2));
            xrfSub=zeros(size(BindingEnergy));
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
                Wsub=reshape(W(index(j,2),index(j,1),:),size(M));
                ExistInd=find(Wsub~=0);
                %% Self-absorption

                BeforeEmit=0;
                I_after=zeros(size(M));
                for tsub=1:length(ExistInd)
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
                        for t_after=1:size(index_after,1)
                            temp_after=temp_after+Lvec_after(t_after)*MU_after{tsub}(index_after(t_after,2),index_after(t_after,1));
                        end %% Attenuation of Flourescent energy emitted from current pixel
                    I_after(tsub)=I_after(tsub)+exp(-temp_after);    
                    end %% End loop for each SSD detector let
                
                end %% End loop for existing fluorescence energy from current pixel

                %% ====================================================================
                RM{index(j,2),index(j,1)}=Lvec(j)*I_incident*(Wsub.*M).*I_after/NumSSDlet;% fluorescence emitted from current pixel
                temp=zeros(size(RMlocal(index(j,2),index(j,1),:)));
                temp(1,1,:)=Lvec(j)*I_incident*(Wsub.*M).*I_after/NumSSDlet;
                RMlocal(index(j,2),index(j,1),:)=RMlocal(index(j,2),index(j,1),:)+temp;
                for tsub=1:length(ExistInd)
                    DetInd=find(Energy(ExistInd(tsub))==BindingEnergy);
                    
                    xrfSub(DetInd)=xrfSub(DetInd)+RM{index(j,2),index(j,1)}(ExistInd(tsub));
                end
            end
            Ltol{n,i}=L;
            p1=GaussianFit1(BindingEnergy,xrfSub);
            if(plotSpec)
                finalfig=figure('name',sprintf('XRF at beam %d with angle %d',i,thetan(n)));
                plot(DetChannel,p1,'r-')
                xlabel('Energy Channel','fontsize',12); ylabel('Intensity','fontsize',12)
                pause;
            end
            XRF{n,i}=p1;
        end
        
    end
end

