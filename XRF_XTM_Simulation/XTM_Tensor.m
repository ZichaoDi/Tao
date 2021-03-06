%%%Simulate XRT of a given object with predifined detector and beam
global plotSpecSingle BeforeEmit N 
global LogScale EmptyBeam synthetic 
plotSpecSingle=0;
more off;
%%-----------------------------------------------
Define_Detector_Beam_Gaussian; %% provide the beam source and Detectorlet
DefineObject_Gaussian; % Produce W, MU_XTM
%%==============================================================
Energy=BindingEnergy;
mtol=prod(m);
eX=ones(m(1),1);
eY=ones(m(2),1);
EmptyBeam=[];
L=sparse(numThetan*(nTau+1),prod(m));
DisR_Simulated=zeros(nTau+1,numThetan);
GlobalInd=cell(numThetan,nTau+1);
fprintf(1,'====== Fluorescence Detector Resolution is %d\n',numChannel);
for n=1:numThetan
    if(mod(n,10)==0)
        fprintf(1,'====== Angle Number  %d of %d: %d\n',n,numThetan,thetan(n));
    end
    theta = thetan(n)/180*pi;
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    DetKnot=DetKnot0*TransMatrix;
    SourceKnot=SourceKnot0*TransMatrix;
    %% =========================================
    Rdis_true=1*mean(I0(:))*ones(nTau+1,1);
    xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
    ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
    plotDisBeam=0;
    if(plotDisBeam)
        finalfig=figure('name',sprintf('XRT with %d beamlets and angle %d',nTau+1,thetan(n)));
    end
    for i=1:nTau+1 %%%%%%%%%========================================================
        % Initialize
        BeforeEmit=1;
        %============================= Plot Grid and Current Light Beam
        if(plotDisBeam)
            % finalfig=figure('name',sprintf('XRF at beam %d with angle %d',i,thetan(n)));
            figure(finalfig)
            hold on;
            subplot(1,2,1)
            [px,py]=meshgrid(linspace(omega(1),omega(2),m(1)+1),linspace(omega(3),omega(4),m(2)+1));
            plot(px,py,'b-',px',py','k-','LineWidth',2);
            axis equal
            dp = -SourceKnot+DetKnot;   
            quiver(SourceKnot(:,1),SourceKnot(:,2),dp(:,1),dp(:,2),0,'color','r');
            drawnow;
        end
        %=================================================================
        %=================================================================
        [index,Lvec,linearInd]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta);
        %%%%%%%%================================================================
        GlobalInd{n,i}=index;
        if(~isempty(index)& norm(Lvec)>0)
            EmptyBeam=[EmptyBeam,(n-1)*numThetan+i];
            currentInd=sub2ind(m,index(:,2),index(:,1));
            L(sub2ind([numThetan,nTau+1],n,i),currentInd)=Lvec;

            Rdis_true(i)=mean(I0(:))*exp(-eX'*(MU_XTM.*reshape(L(sub2ind([numThetan,nTau+1],n,i),:),m))*eY);%%I0*exp(-eX'*(MU_XTM.*reshape(L(n,i,:,:),subm,subn))*eY); %% Discrete case
        end
    end
        DisR_Simulated(:,n)=Rdis_true';
        if(plotSpec)
            finalfig=figure('name',sprintf('XRT with angle %d',thetan(n)));
            subplot(1,2,1)
            plotGrid(xc,omega,[m(2) m(1)]);
            hold on;
            plot(DetKnot(:,1),DetKnot(:,2),'k+-',SourceKnot(:,1),SourceKnot(:,2),'m+-','LineWidth',0.5)
            imagesc((x(1:end-1)+x(2:end))/2,(y(1:end-1)+y(2:end))/2,MU)
            pause(1);
        end
        if(plotDisBeam)
            subplot(1,2,2);
            plot(1:nTau+1,DisR_Simulated(:,n),'r.-')
            xlabel('Beamlet','fontsize',12); ylabel('Intensity','fontsize',12)
        end
end
%%==============================================================
