% Define detector and beam
global m Tol x y omega  dz
global DetChannel numChannel nTau DetKnot0 SourceKnot0 NumSSDlet
Tol=1e-5; %%the threshod to gurantee the beam will cover the whole object
omega=[-2     2    -2     2].*Tol;
m=[current_n current_n]; %Numerical Resolution
%%%%%======================================
DisY=m(1)+1;
DisX=m(2)+1;
x=linspace(omega(1),omega(2),DisX);
y=linspace(omega(3),omega(4),DisY);
dz=[(omega(2)-omega(1))/m(2) (omega(4)-omega(3))/m(1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=atan((omega(4)-omega(3))/(omega(2)-omega(1)));
dTau=(omega(2)-omega(1))/m(1);%%% width of the each discrete beam
Tau=sqrt((omega(2)-omega(1))^2+(omega(4)-omega(3))^2);
nTau=ceil(Tau/dTau)+1;%m(1)+1;% % number of discrete beam%nTau;%
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
%----------------------------------------------------------------
NumSSDlet=5;
SSDlet=[linspace(SSD0(2,1),SSD0(1,1),NumSSDlet)',...
    linspace(SSD0(2,2),SSD0(1,2),NumSSDlet)' ];
