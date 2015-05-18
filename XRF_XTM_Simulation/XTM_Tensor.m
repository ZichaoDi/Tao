%%%Simulate XRT of a given object with predifined detector and beam
global plotSpecSingle BeforeEmit
global LogScale EmptyBeam
plotSpecSingle=0;
more off;
Define_Detector_Beam_Gaussian; %% provide the beam source and Detectorlet
DefineObject_Gaussian; % Produce W, MU_XTM
%%%%%%%==============================================================
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
    
    Rdis=I0*ones(nTau+1,1);
    for i=1:nTau+1 %%%%%%%%%================================================================
        % Initialize
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        BeforeEmit=1;
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
            Rdis(i)=I0*exp(-eX'*(MU_XTM.*reshape(L(n,i,:,:),subm,subn))*eY); %% Discrete case
        end
    end
    DisR(:,n)=Rdis';
end
if(LogScale)
    SigMa_XTM=1./(-log(DisR(:)./I0));
else
    SigMa_XTM=1./DisR(:);
end
SigMa_XTM=ones(size(SigMa_XTM));
