%%=======================================================================
global NumElement  

load PeriodicTable
ScaleM=1e-5;
%%%========================== assign weight matrix for each element in each pixel
Z=[8 30 20 16];% reference sample: Pb La Pd Mo Cu Fe Ca
Z=Z(1:NumElement);
UnitSpectrumSherman_Gaussian; %% Produce BindingEnergy M
%%=======================================================================
NumLines=NumElement;
MU_e=zeros(NumElement,1,1+NumLines);
for i=1: NumLines
    MU_e(i,1,1)=na(i)*CS_TotalBeam(Z(i),1);%calllib('libxrl','CS_Total',Z(i),E0);
    for j=1:NumElement
        MU_e(i,1,j+1)=na(i)*CS_Total(Z(i),1,Z(j));%calllib('libxrl','CS_Total',Z(i),BindingEnergy(j));   only consider K_alpha line
    end
end
%%%%%================== Attenuation Matrix at beam energy


