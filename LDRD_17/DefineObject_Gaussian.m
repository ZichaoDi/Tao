%%=======================================================================
%%% Input: Z: Atomic Numbers of Existing Elements in the object
%%%        E0: Beam Energy
%%%        m: Resolution of object
%%% Output: W: Concentration fraction of each Z in each pixel with
%%%            dimension m(1) by m(2) by cardinality(Z)
%%%         O: Object specified by W.*Z
%%%         MU: Attenuation matrix of O
%%=======================================================================
global x y dz m omega MU_e Z 
global NumElement Element 
global slice 

load PeriodicTable
%%%%%======================================
x=linspace(omega(1),omega(2),m(1)+1);
y=linspace(omega(3),omega(4),m(2)+1);

%%%=========== assign weight matrix for each element in each pixel
%---------------------------
UnitSpectrumSherman_real; %% Produce BindingEnergy M
clear Line ElementDensity LineEnergy CS_FluoLine CS_Total CS_TotalBeam
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,1,NumElement);
MU_XTM=sum(1e0*W.*repmat(MUe,[m(1),m(2),1]),3);
%%%%% ====================================================================

%%=================== Picture the object based on atomic number and attenuation

