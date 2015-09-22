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
global x y m omega dz AbsorbScale MU_e Z
global XTMscale NumLines NumElement
global MUe N synthetic 

load PeriodicTable
AbsorbScale=1e-0;
ScaleM=1.5e-4;
%%%%%======================================
DisY=m(1)+1;
DisX=m(2)+1;
x=linspace(omega(1),omega(2),DisX);
y=linspace(omega(3),omega(4),DisY);
dz=[(omega(2)-omega(1))/m(2) (omega(4)-omega(3))/m(1)];
[X,Y] = meshgrid((x(1:end-1)+x(2:end))./2,(y(1:end-1)+y(2:end))./2);
center=[0 0];
%%%========================== the grids of object
xc = getNodalGrid(omega,[m(2) m(1)]);
%%%========================== assign weight matrix for each element in each pixel
if(synthetic)
Z = [6 8 14 20 26];% Golosio's Sample
else
% Z=[14 29 30 74 79];% Glass Rod
Z=[14 16 17 19 20 22 23 24 25 26 28 29 30 31 80 33 34 35 92 37 38 39 40];% Complete Seed
end
%---------------------------
NumElement=length(Z);
W=ones(N(level),N(level),NumElement);
if(level==1)
    if(synthetic)
    CreateCircle; %% Golosio's sample 
    else
         for tsub=1:NumElement
          W(:,:,tsub)=abs(flipud(permute(iR_num(:,:,slice(tsub)),[2 1 3])));%tsub*2e-1;% Seed
        % % W(:,:,tsub)=abs(fliplr(rot90(permute(iR(tsub,:,:),[2 3 1]))));%tsub*2e-1;% GlassRod
         end
    end
end
UnitSpectrumSherman_Gaussian; %% Produce BindingEnergy M
%%=======================================================================
NumLines=NumElement;
MU_e=zeros(NumElement,1,1+NumLines);
for i=1: NumLines
    MU_e(i,1,1)=calllib('libxrl','CS_Total',Z(i),E0);% na(i)*CS_TotalBeam(Z(i),1);%
    for j=1:NumElement
        MU_e(i,1,j+1)=na(i)*calllib('libxrl','CS_Total',Z(i),BindingEnergy(j));%CS_Total(Z(i),1,Z(j));%   only consider K_alpha line
    end
end
%%%%%================== Attenuation Matrix at beam energy
MU_e=MU_e.*AbsorbScale; %% Discrete Scale
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,1,NumElement);
MU=sum(W.*repmat(MUe,[m(1),m(2),1]),3);
%%=======smooth data
% MU=BlurGaussian(MU);
%{
yy=reshape(smooth(MU(:)),m(1),m(2));
MU=reshape(smooth(reshape(yy',prod(m),1)),m(1),m(2));
%}
%%===================
XTMscale=1e0;
MU_XTM=MU.*XTMscale;
% MU_XTM=ir';

%%%%% =================== Attenuation Matrix at flourescence energy (Corrected Attenuation)
MU_after=cell(NumElement,1);
for i=1:NumElement
    MU_after{i}=sum(W.*repmat(reshape(MU_e(:,1,i+1),1,1,NumElement),[m(1),m(2),1]),3);
end
%%%%% ====================================================================
%%=================== Picture the object based on atomic number and attenuation
if(PlotObject)
    figureObject(W,Z,m,NumElement,MU_e,0)
end

%%=======================================================================
