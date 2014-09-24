% Define detector and beam 
global m ;
 global DetChannel numChannel nTau DetKnot0 SourceKnot0 NumSSDlet
omega=[-2     2    -2     2];
  m=[50 50];
 alpha=atan((omega(4)-omega(3))/(omega(2)-omega(1)));
dTau=(omega(2)-omega(1))/m(2);%0.5;%%% width of the each discrete beam
Tau=sqrt((omega(2)-omega(1))^2+(omega(4)-omega(3))^2);
nTau=ceil(Tau/dTau)+2; % number of discrete beam
tol=10; %%the threshod to gurantee the beam will cover the whole object
detS0=[Tau/2*tan(alpha)+tol, ceil(Tau/2)]; %initiate transmission detector location
detE0=[Tau/2*tan(alpha)+tol,-ceil(Tau/2)];
knot=linspace(detS0(2),detE0(2),nTau+1)';
DetKnot0=[repmat(detS0(1),size(knot)),knot];%% transmission detectorlet knot points
%%%=====================================================
SourceS0=[-Tau/2*tan(alpha)-tol, ceil(Tau/2)];%initiate beam source
SourceE0=[-Tau/2*tan(alpha)-tol,-ceil(Tau/2)];
knot=linspace(SourceS0(2),SourceE0(2),nTau+1)';
SourceKnot0=[repmat(SourceS0(1),size(knot)),knot];%% source knot points
%%%======================Define Energy Channel of the fluorescence detector
SSD0=[detE0-[0,tol]; SourceE0-[0,tol]];
numChannel=5;
DetScaleXRF=5;
DetChannel=linspace(0,DetScaleXRF,numChannel)';
NumSSDlet=2;
SSDlet=[linspace(SSD0(2,1),SSD0(1,1),NumSSDlet)',...
            linspace(SSD0(2,2),SSD0(1,2),NumSSDlet)' ];
