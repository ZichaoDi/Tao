% Define detector and beam 
 global m DetChannel numChannel nTau DetKnot0 SourceKnot0 dz

dz=Dis/m(1); %%% width of the each discrete beam
alpha=atan(m(2)/m(1));
dTau=dz;
Tau=sqrt(m(1)^2+m(2)^2);
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
SSD0=[detE0; SourceE0];
numChannel=10;
DetChannel=linspace(0,100,numChannel)';
NumSSDlet=10;
SSDlet=[linspace(SSD0(2,1),SSD0(1,1),NumSSDlet)',...
            linspace(SSD0(2,2),SSD0(1,2),NumSSDlet)' ];
