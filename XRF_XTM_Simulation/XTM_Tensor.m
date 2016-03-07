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
if(level==1)
    DisR=zeros(nTau+1,numThetan);
    EmptyBeam=[];
end
if(level==1)
   L=sparse(numThetan*(nTau+1),prod(m));
else
   L=sparse(numThetan*(nTau+1),prod(m));
   % L_H=sparse(numThetan*(nTau+1),prod(m));
end
GlobalInd=cell(numThetan,nTau+1);
fprintf(1,'====== Fluorescence Detector Resolution is %d\n',numChannel);
pert_drift=0;
pert_angle=0;
drift_angle=pert_angle*0.2*rand(size(thetan));
for n=1:numThetan
    theta=thetan(n)/180*pi+drift_angle(n);
    if(mod(n,10)==0)
        fprintf(1,'====== Angle Number  %d of %d: %d\n',n,numThetan,thetan(n));
    end
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    driftx=pert_drift*1e-3*mean(DetKnot0(:,1))*rand(size(DetKnot0(:,1)));
    drifty=pert_drift*1e-3*mean(DetKnot0(:,2))*rand(size(DetKnot0(:,2)));
    DetKnot=DetKnot0*TransMatrix+[driftx,drifty];
    SourceKnot=SourceKnot0*TransMatrix+[driftx,drifty];
    Rdis=1*I0*ones(nTau+1,1);
    xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
    ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
    for i=1:nTau+1 %%%%%%%%%================================================================
        %if(level==1)
        % Initialize
        BeforeEmit=1;
        %============================= Plot Grid and Current Light Beam
        plotDisBeam=0;
        if(plotDisBeam)
            finalfig=figure('name',sprintf('XRF at beam %d with angle %d',i,thetan(n)));
            subplot(1,2,1)
            plotGrid(xc,omega,[m(2) m(1)]);
            hold on;
            plot(DetKnot(:,1),DetKnot(:,2),'k+-',SourceKnot(:,1),SourceKnot(:,2),'m+-','LineWidth',0.5)
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
        %=================================================================
        [index,Lvec,linearInd]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta);
        %%%%%%%%================================================================
        GlobalInd{n,i}=index;
        if(~isempty(index)& norm(Lvec)>0)
            EmptyBeam=[EmptyBeam,(n-1)*numThetan+i];
            currentInd=sub2ind(m,index(:,2),index(:,1));
            L(sub2ind([numThetan,nTau+1],n,i),currentInd)=Lvec;
            Rdis(i)=I0*exp(-eX'*(MU_XTM.*reshape(L(sub2ind([numThetan,nTau+1],n,i),:),m))*eY);%%I0*exp(-eX'*(MU_XTM.*reshape(L(n,i,:,:),subm,subn))*eY); %% Discrete case
        end

        % else
        %     L_H(sub2ind([numThetan,nTau+1],n,i),:)= ...
        %     L(sub2ind([numThetan,nTau+1],n,i),:)*IhH{level}';
        % end
    end
        if(plotSpec)
            finalfig=figure('name',sprintf('XRT with angle %d',thetan(n)));
            subplot(1,2,1)
            plotGrid(xc,omega,[m(2) m(1)]);
            hold on;
            plot(DetKnot(:,1),DetKnot(:,2),'k+-',SourceKnot(:,1),SourceKnot(:,2),'m+-','LineWidth',0.5)
            imagesc((x(1:end-1)+x(2:end))/2,(y(1:end-1)+y(2:end))/2,MU)
            subplot(1,2,2);
            plot(1:nTau+1,Rdis,'r-')
            xlabel('Detector Channel','fontsize',12); ylabel('Intensity','fontsize',12)
            pause(1);
        end
     if(level==1)
         DisR(:,n)=Rdis';
     end
end
if(level>1)
L=L;%L_H;
end
%%==============================================================
if(~synthetic)
    DisR_Simulated=DisR;
    DisR=squeeze(data_xrt)';% reshape(dd,size(dd,1)/numThetan,numThetan);%
end
if(LogScale)
    SigMa_XTM=1./(-log(DisR(:)));
else
    SigMa_XTM=1./DisR(:);
end
SigMa_XTM=ones(size(SigMa_XTM));
