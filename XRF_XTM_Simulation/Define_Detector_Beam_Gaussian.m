% Define detector and beam 
global m Tol thetan
global DetChannel numChannel nTau DetKnot0 SourceKnot0 NumSSDlet 
Tol=1e-5; 
omega=[-2     2    -2     2].*Tol;
m=[current_n current_n]; %Numerical Resolution
alpha=atan((omega(4)-omega(3))/(omega(2)-omega(1)));
dTau=(omega(2)-omega(1))/N(1);%%% width of each discrete beam
Tau=sqrt((omega(2)-omega(1))^2+(omega(4)-omega(3))^2);
if(synthetic)
    nTau=m(1)+1;%ceil(Tau/dTau)+1;% % number of discrete beam%nTau;%
else
    nTau=size(data_h,3)-1;
end
tol1=eps^(1/2); % the threshod to gurantee the beam will cover the whole object
%=============initiate transmission detector location
detS0=[Tau/2*tan(alpha)+tol1*Tol, Tau/2+tol1*Tol]; 
detE0=[Tau/2*tan(alpha)+tol1*Tol,-Tau/2-tol1*Tol];
knot=linspace(detS0(2),detE0(2),nTau+1)';
DetKnot0=[repmat(detS0(1),size(knot)),knot];%% transmission detectorlet knot points
SourceS0=[-Tau/2*tan(alpha)-tol1*Tol, Tau/2+tol1*Tol];%initiate beam source
SourceE0=[-Tau/2*tan(alpha)-tol1*Tol,-Tau/2-tol1*Tol];

knot=linspace(SourceS0(2),SourceE0(2),nTau+1)';
SourceKnot0=[repmat(SourceS0(1),size(knot)),knot];%% source knot points
SSD0=[detE0-[0,Tol]; SourceE0-[0,Tol]];
%%%=============Define Energy Channel of the fluorescence detector
if(synthetic)
    DetScaleXRF=10;
    DetChannel=linspace(0,DetScaleXRF,numChannel)';
else
    load(['DetChannel_',sample,'.mat']); %-- APS real fluorescence detector energy channel
    numChannel=length(DetChannel);
    if(truncChannel)
        truncInd=find(DetChannel> 5 & DetChannel <9);
        DetChannel=DetChannel(truncInd);
        numChannel=length(truncInd);
    end
    DetChannel_raw=DetChannel;
    if (DecomposedElement)
        numChannel=length(slice);
        DetChannel=[1:numChannel]';
    end
end
%%%=========== Define number of flying paths of fluorescence photons
NumSSDlet=5;
SSDlet=[linspace(SSD0(2,1),SSD0(1,1),NumSSDlet)',...
            linspace(SSD0(2,2),SSD0(1,2),NumSSDlet)' ];
%%%=========== Assign Projection Angles;
thetan=linspace(363,abs(183*(angleScale)-363),numThetan);% must be positive.
if(strcmp(sample,'Rod'))
    thetan=linspace(-180,180,numThetan)+360;% must be positive.
end
