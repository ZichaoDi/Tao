% Net intensity = I0*c_e*omega_e*(1-1/JumpFactor)*mu*T
% c_e: element concentration , for single element, c_e=1
% omega_e: fluorescence yield
% mu: linear attenuation coefficient at I0
% T: sample thickness, for single element, assume constant T=1 for all
% n_a: the number density of atoms
% mu=n_a*CrossSection
NA=6.02e23;%Avogadro's number
loadlibrary('/opt/local/lib/libxrl.dylib','/opt/local/include/xraylib/xraylib.h');
load PeriodicTable
Line=2; %% Transition Line
shell=0;  %% Shell type
BindingEnergy=zeros(NumElement*length(Line),1);
M=zeros(NumElement*length(Line),1);
i=1;
T=1;
if(plotUnit)
    figure('name', 'Unit Spectrum');
end
fprintf(1,'Element     AtomicNumber    Line     FluorescenceEnergy     Intensity\n');
while(i<=NumElement)
    density=calllib('libxrl','ElementDensity',Z(i));
    A=calllib('libxrl','AtomicWeight',Z(i));
    na(i)=density*NA/A*Z(i)*1e-24;
    j=1;
    while (j<=length(Line))
        new_energy=calllib('libxrl','LineEnergy',Z(i),-Line(j));
        mu=na(i)*calllib('libxrl','CS_Total',Z(i),E0);
        Jump=calllib('libxrl','JumpFactor',Z(i),shell);
        w=calllib('libxrl','FluorYield',Z(i),Line(j)-1);
        c=1; %%% weight fraction, 1: single element, ~1: compound
        intensity=E0*c*w*(1-1/Jump)*mu*T;
        fprintf('    %s           %d          %i          %5.2f               %5.2f\n', Element{Z(i)}, Z(i),-Line(j), new_energy, intensity);
        BindingEnergy(length(Line)*(i-1)+j)=new_energy;
        M(length(Line)*(i-1)+j)=intensity;
        if(intensity~=0)
            theta=intensity^(1/3);
            miu=new_energy;
            p=intensity/(theta*sqrt(2*pi))*exp(-(DetChannel-miu).^2./(2*theta^2));
            if(plotUnit)
            plot(DetChannel,p,'r-')
            title('Gaussian spectrum for each element')
            hold on;
            end
        end
        %%%%===================================================================
        
        j=j+1;
    end
    fprintf('=====================================================\n')
    i=i+1;
end

