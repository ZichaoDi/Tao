% Net intensity = I0*c_e*omega_e*(1-1/JumpFactor)*mu*T
% c_e: element concentration , for single element, c_e=1
% omega_e: fluorescence yield
% mu: linear attenuation coefficient at I0
% T: sample thickness, for single element, assume constant T=1 for all
% n_a: the number density of atoms
% mu=n_a*CrossSection
global TakeLog
E0=40;
E2I=1/3e8/6.62e-34*1e-15;
I0=E0*E2I;
Z=[42 29 26 20];%[82 57 46 42 29 26 20];% reference sample: Pb La Pd Mo Cu Fe Ca
NumElement=length(Z);
NA=6.02e23;%Avogadro's number
if(ismac)
    loadlibrary('/opt/local/lib/libxrl.dylib','/opt/local/include/xraylib/xraylib.h');
else
    loadlibrary('/homes/wendydi/Documents/lib/libxrl.so','/homes/wendydi/Documents/include/xraylib/xraylib.h');
end
load PeriodicTable
Line=[-0];% -1 -2 -3]; %% Transition Line
shell=0;  %% Shell type
BindingEnergy=zeros(NumElement*length(Line),1);
M=zeros(NumElement*length(Line),numChannel);

i=1;
T=1;
if(plotUnit)
    figure('name', 'Unit Spectrum');
end
fprintf(1,'Element     AtomicNumber    Line     FluorescenceEnergy     Intensity\n');
%%==set up Gaussian Distribution
mu=0;
sigma=(DetChannel(2)-DetChannel(1))/2.35;
G=fft(1*exp(-(DetChannel-mu).^2./(2*sigma^2))/(sigma*sqrt(2*pi)));
%%======================================================================
while(i<=NumElement)
    PurePeak=0.*DetChannel;
    density=calllib('libxrl','ElementDensity',Z(i));
    A=calllib('libxrl','AtomicWeight',Z(i));
    na(i)=density*NA/A*Z(i)*1e-24;
    j=1;
    while (j<=length(Line))
        new_energy=calllib('libxrl','LineEnergy',Z(i),-Line(j));
        %         mu=na(i)*calllib('libxrl','CS_Total',Z(i),E0);
        %         Jump=calllib('libxrl','JumpFactor',Z(i),shell);
        %         w=calllib('libxrl','FluorYield',Z(i),Line(j)-1);
        %         c=1; %%% weight fraction, 1: single element, ~1: compound
        %         intensity=I0*c*w*(1-1/Jump)*mu*T;
        intensity=I0*calllib('libxrl','CS_FluorLine',Z(i),-Line(j),E0);
        fprintf('    %s           %d          %i          %5.2f               %5.2f\n', Element{Z(i)}, Z(i),-Line(j), new_energy, intensity);
        BindingEnergy(length(Line)*(i-1)+j)=new_energy;
        [~,adj]=min(abs(repmat(new_energy,size(DetChannel))-DetChannel));
        PurePeak(adj)=PurePeak(adj)+intensity;
        
        j=j+1;
    end
    %%%%=============================== Start Gaussian Convolution
    M(i,:)=ifft(fft(PurePeak).*G);
    if(plotUnit)
        cmap=colormap(lines);
        cmap = cmap(1:NumElement,:);
        plot(DetChannel,M(i,:),'color',cmap(i,1:3),'LineStyle',':','LineWidth',1.5);
        Legend{i}=sprintf('%s',Element{Z(i)});
        legend(Legend)
        title('Gaussian spectrum for each element')
        hold on;
    end
    
    %%%%===================================================================
    fprintf('=====================================================\n')
    i=i+1;
end
M=M.*1e-9;
if(TakeLog)
M=abs(M);%1e-300;
end



