% Define detector and beam 
global m Tol thetan
global DetChannel numChannel nTau DetKnot0 SourceKnot0 NumSSDlet 
% Tol=1e-5; 
% omega=[-2     2    -2     2].*Tol;
m=[current_n current_n]; %Numerical Resolution
alpha=atan((omega(4)-omega(3))/(omega(2)-omega(1)));
dTau=(omega(2)-omega(1))/N(1);%%% width of each discrete beam
Tau=sqrt((omega(2)-omega(1))^2+(omega(4)-omega(3))^2)-Tol;% omega(2)-omega(1);%
if(synthetic)
    nTau=m(1)+1;%ceil(Tau/dTau)+1;% % number of discrete beam%nTau;%
else
    nTau=size(data_h,3)-1;
end
tol1=0;%eps^(1/2); % the threshod to gurantee the beam will cover the whole object
%=============initiate transmission detector location
detS0=[Tau/2*tan(alpha)+tol1*Tol, Tau/2+tol1*Tol]; 
detE0=[Tau/2*tan(alpha)+tol1*Tol,-Tau/2-tol1*Tol];
knot=linspace(detS0(2),detE0(2),nTau+1)';
DetKnot0=[repmat(detS0(1),size(knot)),knot];%% transmission detectorlet knot points
SourceS0=[-Tau/2*tan(alpha)-tol1*Tol, Tau/2+tol1*Tol];%initiate beam source
SourceE0=[-Tau/2*tan(alpha)-tol1*Tol,-Tau/2-tol1*Tol];

knot=linspace(SourceS0(2),SourceE0(2),nTau+1)';
SourceKnot0=[repmat(SourceS0(1),size(knot)),knot];%% source knot points
DetKnot0=DetKnot0(end:-1:1,:);
SourceKnot0=SourceKnot0(end:-1:1,:);
SSD0=[detE0-[0,Tol]; SourceE0-[0,Tol]];
%%%=============Define Energy Channel of the fluorescence detector
if(synthetic)
    DetScaleXRF=10;
    DetChannel=linspace(0,DetScaleXRF,numChannel)';
else
    % ==== APS real fluorescence detector energy channel
    load(['DetChannel_',sample,'.mat']); 
    numChannel_raw=length(DetChannel);
    DetChannel_raw=DetChannel;
    % ==== Decomposed channel
    numChannel_decom=size(iR_num,3);
    DetChannel_decom=[1:numChannel_decom]';
end
%%%=========== Define number of flying paths of fluorescence photons
NumSSDlet=5;
SSDlet=[linspace(SSD0(2,1),SSD0(1,1),NumSSDlet)',...
            linspace(SSD0(2,2),SSD0(1,2),NumSSDlet)' ];
%%%=========== Assign Projection Angles;
thetan=linspace(363,abs(183*(angleScale)-363),numThetan);% must be positive.
if(strcmp(sample,'Rod'))
    thetan=thetan_real;%linspace(-180,180,numThetan)+360;% must be positive.
end
