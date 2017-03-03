% Define detector and beam 
global m Tol dTau thetan
global DetChannel numChannel nTau DetKnot0 SourceKnot0 NumSSDlet 
alpha=atan((omega(4)-omega(3))/(omega(2)-omega(1)));
Tau= omega(2)-omega(1);
if(synthetic)
    nTau=1*ceil(sqrt(2*prod(m)));%m(1)-1;% % number of discrete beam%nTau;%
    tol1=1/2*N(1);
else
    tol1=0;
end
%=============initiate transmission detector location
detS0=[Tau/2*tan(alpha)+tol1*dz(1), Tau/2+tol1*dz(1)]; 
detE0=[Tau/2*tan(alpha)+tol1*dz(1),-Tau/2-tol1*dz(1)+delta_d0];
dTau=abs(-Tau-2*tol1*dz(1)+delta_d0)/(nTau+1);%%% width of each discrete beam
knot=linspace(detS0(2),detE0(2),nTau+1)';
DetKnot0=[repmat(detS0(1),size(knot)),knot];%% transmission detectorlet knot points
SourceS0=[-Tau/2*tan(alpha)-tol1*dz(1), Tau/2+tol1*dz(1)];%initiate beam source
SourceE0=[-Tau/2*tan(alpha)-tol1*dz(1),-Tau/2-tol1*dz(1)+delta_d0];

knot=linspace(SourceS0(2),SourceE0(2),nTau+1)';
SourceKnot0=[repmat(SourceS0(1),size(knot)),knot];%% source knot points
DetKnot0=DetKnot0(end:-1:1,:);
SourceKnot0=SourceKnot0(end:-1:1,:);
SSD0=[detE0-[0,Tol]; SourceE0-[0,Tol]];
%%%=============Define Energy Channel of the fluorescence detector
% ==== APS real fluorescence detector energy channel
load(['DetChannel_','Rod','.mat']); 
numChannel=length(DetChannel);
DetChannel_raw=DetChannel;
numChannel_raw=length(DetChannel);
% ==== Decomposed channel
numChannel_decom=NumElement;
DetChannel_decom=[1:numChannel_decom]';
%%%=========== Define number of flying paths of fluorescence photons
NumSSDlet=5;
SSDlet=[linspace(SSD0(2,1),SSD0(1,1),NumSSDlet)',...
            linspace(SSD0(2,2),SSD0(1,2),NumSSDlet)' ];
%%%=========== Assign Projection Angles;
thetan=linspace(0,360,numThetan);%linspace(363,abs(183*(angleScale)-363),numThetan);% must be positive.
if(numThetan==1)
    thetan=90;
end
if(strcmp(sample,'Rod'))
    thetan=thetan_real;%linspace(-180,180,numThetan)+360;% must be positive.
end