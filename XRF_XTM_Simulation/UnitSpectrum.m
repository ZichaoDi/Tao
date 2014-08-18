%% Calculate the unit spectrum of existing element with beam energy E0 without Gaussian fit
%% Input: Z, E0
%% Output: BindingEnergy: Binding energy needed to excite Z
%%         M: Unit intensity of Z corresponding to BindingEnergy

load PeriodicTable
Line=2; %% Transition Line: 2-Kalpha2
load([num2str(E0),'ElementStrKalpha',num2str(Line),'.mat']);
flux=1e5; %% Incident photon flux
G=1e-4; % geometry factor determined by the characteristics of the detector and its position
        % with respect to the beam
i=1;
new_energy=E0;
thickness=0.01; %%Thickness of element
BindingEnergy=zeros(NumElement*length(Line),1);
M=zeros(NumElement*length(Line),1);
fprintf(1,'Element     AtomicNumber    Line     FluorescenceEnergy     Intensity\n');
while(i<=length(Z))
    Density=density(Z(i));

    mu_0=mu0(Z(i)); %% Attenuation Coefficient at beam energy
    j=1;
    while (j<=length(Line))
        new_energy=Kalpha2Energy(Z(i)); %% Binding Eenergy of Z
        mu_1=mu1(Z(i)); %% Attenuation Coefficient at Fluorescene Energy
        Q=CS(Z(i)); %% Cross Section
        w=1;%weight fraction of compound
        chi=mu_0+mu_1;
        A_corr=1-exp(-chi*Density*thickness);
        intensity=flux*G*Q*w*A_corr;
        
        if(intensity~=0)
            fprintf('    %s           %d          %i          %5.2f               %5.2f\n', Element{Z(i)}, Z(i),-Line(j), new_energy, intensity);
        end
        BindingEnergy(length(Line)*(i-1)+j)=new_energy;
        M(length(Line)*(i-1)+j)=intensity;
        
        %%%%===================================================================
        j=j+1;
    end
    fprintf('=====================================================\n')
    i=i+1;
end
