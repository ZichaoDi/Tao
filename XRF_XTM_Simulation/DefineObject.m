% function [O,MU,W]=DefineObject(Z,E0)
%%=======================================================================
%%% Input: Z: Atomic Numbers of Existing Elements in the object 
%%%        E0: Beam Energy
%%%        m: Resolution of object
%%% Output: W: Concentration fraction of each Z in each pixel with 
%%%            dimension m(1) by m(2) by cardinality(Z)
%%%         O: Object specified by W.*Z
%%%         MU: Attenuation matrix of O
%%=======================================================================
global x y m omega NumElement
Z=[30 20 35 26];
NumElement=length(Z);
E0=20;
UnitSpectrumSherman;
load([num2str(E0),'ElementStrKalpha2','.mat']);
m=[3 3]; 
%%%%%======================================
Dis=m(1)+1;
omega=[-Dis/2 Dis/2 -Dis/2 Dis/2];
x=linspace(omega(1),omega(2),Dis);
y=linspace(omega(3),omega(4),Dis);
[X,Y] = meshgrid((x(1:end-1)+x(2:end))./2,(y(1:end-1)+y(2:end))./2);
center=[0 0];
% MU=(X-center(1)).^2+(Y-center(2)).^2;
%%%========================== plot the grids of object
if(plotTravel)
xc = getNodalGrid(omega,m);
fig=figure('name','grids & beam line');
plotGrid(xc,omega,m); hold on;
end
%%%%%================== Attenuation Matrix at beam energy
MU=zeros(m);
MU(1,1)=mu0(Z(1));
MU(1,2)=mu0(Z(2));
MU(1,3)=mu0(Z(1));
MU(2,1)=mu0(Z(2));
MU(2,2)=0.5*(mu0(Z(3))+mu0(Z(4)));
MU(2,3)=mu0(Z(2));
MU(3,1)=mu0(Z(3));
MU(3,2)=mu0(Z(1));
MU(3,3)=mu0(Z(2));
MU=MU.*1e-1;
%%%%% =================== Attenuation Matrix at flourescence energy (Corrected Attenuation)
MU_after=cell(NumElement,1);
for i=1:NumElement
MU_after{i}(1,1)=na(1)*calllib('libxrl','CS_Total',Z(1),BindingEnergy(i));
MU_after{i}(1,2)=na(2)*calllib('libxrl','CS_Total',Z(2),BindingEnergy(i));
MU_after{i}(1,3)=na(1)*calllib('libxrl','CS_Total',Z(1),BindingEnergy(i));
MU_after{i}(2,1)=na(2)*calllib('libxrl','CS_Total',Z(2),BindingEnergy(i));
MU_after{i}(2,2)=0.5*(na(3)*calllib('libxrl','CS_Total',Z(3),BindingEnergy(i))+na(4)*calllib('libxrl','CS_Total',Z(4),BindingEnergy(i)));
MU_after{i}(2,3)=na(2)*calllib('libxrl','CS_Total',Z(2),BindingEnergy(i));
MU_after{i}(3,1)=na(2)*calllib('libxrl','CS_Total',Z(3),BindingEnergy(i));
MU_after{i}(3,2)=na(1)*calllib('libxrl','CS_Total',Z(1),BindingEnergy(i));
MU_after{i}(3,3)=na(2)*calllib('libxrl','CS_Total',Z(2),BindingEnergy(i));
MU_after{i}=MU_after{i}.*1e-1;
end

W=zeros(m(1),m(2),NumElement); 
O=zeros(m(1),m(2),NumElement); 

W(1,1,:)=[1 0 0 0];
W(1,2,:)=[0 1 0 0];
W(1,3,:)=[1 0 0 0];
W(2,1,:)=[0 1 0 0];
W(2,2,:)=   [0 0 0.5 0.5] ;
W(2,3,:)=  [0 1 0 0];
W(3,1,:)=   [0 0 1 0] ;
W(3,2,:)=   [1 0 0 0] ;
W(3,3,:)=  [0 1 0 0]; 

for e=1:NumElement
    O(:,:,e)=W(:,:,e).*Z(e);
end
% figure('name','AttenuationMatrix');
% imagesc(MU); axis xy
figure('name','Object');
for e=1:NumElement
    subplot(2,2,e);
    image(O(:,:,e));
    axis xy
end