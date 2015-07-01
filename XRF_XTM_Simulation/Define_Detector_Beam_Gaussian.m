% Define detector and beam 
global m Tol thetan nTau_level N
global  DecomposedElement DetChannel numChannel nTau DetKnot0 SourceKnot0 NumSSDlet 
Tol=1e-4; %%the threshod to gurantee the beam will cover the whole object
omega=[-2     2    -2     2].*Tol;
% Tol=1e-2;
% omega=[min(x) max(x) min(x) max(x)]*Tol;
m=[current_n current_n]; %Numerical Resolution
% subm=1;
% m=subm.*[3 3];
alpha=atan((omega(4)-omega(3))/(omega(2)-omega(1)));
dTau=(omega(2)-omega(1))/(N(1));%%% width of the each discrete beam
Tau=sqrt((omega(2)-omega(1))^2+(omega(4)-omega(3))^2);
% if(current_n==N(1))
nTau=size(data,1)-1;%ceil(Tau/dTau)+1;%m(1)+1;% % number of discrete beam%nTau;%
% else
% nTau=nTau_level(current_n==N);
% end
tol1=eps^(1/2);
detS0=[Tau/2*tan(alpha)+tol1*Tol, Tau/2+tol1*Tol]; %initiate transmission detector location
detE0=[Tau/2*tan(alpha)+tol1*Tol,-Tau/2-tol1*Tol];
knot=linspace(detS0(2),detE0(2),nTau+1)';
DetKnot0=[repmat(detS0(1),size(knot)),knot];%% transmission detectorlet knot points
%%%===================================================0==
SourceS0=[-Tau/2*tan(alpha)-tol1*Tol, Tau/2+tol1*Tol];%initiate beam source
SourceE0=[-Tau/2*tan(alpha)-tol1*Tol,-Tau/2-tol1*Tol];

knot=linspace(SourceS0(2),SourceE0(2),nTau+1)';
SourceKnot0=[repmat(SourceS0(1),size(knot)),knot];%% source knot points
%%%======================Define Energy Channel of the fluorescence detector
SSD0=[detE0-[0,Tol]; SourceE0-[0,Tol]];
%--------------------------------- APS real fluorescence detector energy channel

load DetChannel
numChannel=length(DetChannel);
if (DecomposedElement)
numChannel=length(slice);
DetChannel=[0:numChannel]';
end
%----------------------------------------------------------------
% numChannel=2;
% DetScaleXRF=numChannel;
% DetChannel=linspace(0,DetScaleXRF,numChannel)';

%----------------------------------------------------------------
NumSSDlet=5;
SSDlet=[linspace(SSD0(2,1),SSD0(1,1),NumSSDlet)',...
            linspace(SSD0(2,2),SSD0(1,2),NumSSDlet)' ];
% Acquire2Daps;
thetan=[0 45];%linspace(0,2*180,numThetan);
% thetan=thetan-90;
thetan=mod(thetan+360,360);%[0 60];%[1 60];%[1:40:180];% Projection Angles, has to be positive.
% subTheta=1:length(thetan);
% thetan=thetan(subTheta);
