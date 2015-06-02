%%%Simulate XRT of a given object with predifined detector and beam
global plotSpecSingle BeforeEmit
global LogScale EmptyBeam
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
DisR=zeros(nTau+1,numThetan);
Ltol=cell(numThetan,nTau+1);
L=zeros(numThetan,nTau+1,m(1),m(2));
GlobalInd=cell(numThetan,nTau+1);
fprintf(1,'====== Fluorescence Detector Resolution is %d\n',numChannel);
EmptyBeam=[];
for n=1:numThetan
    theta=thetan(n)/180*pi;
    fprintf(1,'====== Angle Number  %d of %d: %d\n',n,numThetan,thetan(n));
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    DetKnot=DetKnot0*TransMatrix;
    SourceKnot=SourceKnot0*TransMatrix;
    Rdis=0*I0*ones(nTau+1,1);
    for i=1:nTau+1 %%%%%%%%%================================================================
        % Initialize
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        BeforeEmit=1;
        %============================= Plot Grid and Current Light Beam
        if(plotSpec & i==1)
            finalfig=figure('name',sprintf('XRF at beam %d with angle %d',i,thetan(n)));
            subplot(1,2,1);
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
        
        if(~isempty(index))
            EmptyBeam=[EmptyBeam,(n-1)*numThetan+i];
            for j=1:size(index,1)
                CurrentCellCenter=[(index(j,1)-1/2)*dz(1)-abs(omega(1)),(index(j,2)-1/2)*dz(2)-abs(omega(3))];
                currentInd=sub2ind(m,index(j,2),index(j,1));
                L(n,i,index(j,2),index(j,1))=Lvec(j);
                
            end
            [~,~,subm,subn]=size(L(n,i,:,:));
            Rdis(i)=eX'*(MU_XTM.*reshape(L(n,i,:,:),subm,subn))*eY;%%I0*exp(-eX'*(MU_XTM.*reshape(L(n,i,:,:),subm,subn))*eY); %% Discrete case
        end
    end
    DisR(:,n)=Rdis';
end
%  DisR_real=reshape(XRF(:,slice,5,:),nTau+1,numThetan);
%%==============================================================
DisR_real=DisR;
%  DisR=data;%DisR_real;
if(LogScale)
    SigMa_XTM=1./(-log(DisR(:)./I0));
else
    SigMa_XTM=1./DisR(:);
end
SigMa_XTM=ones(size(SigMa_XTM));
