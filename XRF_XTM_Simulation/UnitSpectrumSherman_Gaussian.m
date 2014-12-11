% Net intensity = I0*c_e*omega_e*(1-1/JumpFactor)*mu*T
% c_e: element concentration , for single element, c_e=1
% omega_e: fluorescence yield
% mu: linear attenuation coefficient at I0
% T: sample thickness, for single element, assume constant T=1 for all
% n_a: the number density of atoms
% mu=n_a*CrossSection
% K_shell:0 L1:1 L2:2 L3:3 M1:4 M2:5 M3:6 M4:7 M5:8
% KA_LINE:0 KB:1 LA:2 LB:3
global TakeLog
E0=40;
E2I=1/3e8/6.62e-34*1e-15;
I0=E0*E2I;
Z=20:70;%];%[29 30 74 79];%% 42 29 26 ];%20 49 57 46];% reference sample: Pb La Pd Mo Cu Fe Ca
% Z=Z(1:NumElement);
NumElement=length(Z);
NA=6.02e23;%Avogadro's number
% if(ismac)
%     loadlibrary('/opt/local/lib/libxrl.dylib','/opt/local/include/xraylib/xraylib.h');
% else
%     loadlibrary('/homes/wendydi/Documents/lib/libxrl.so','/homes/wendydi/Documents/include/xraylib/xraylib.h');
% end
load AtomicWeight
load(['xRayLib',num2str(E0),'.mat'])
load PeriodicTable
Line=1:36;%[-0 -1 -2 -3]; %% Transition Line, detailed defination see xraylib-lines.h
shell=0;  %% Shell type
BindingEnergy=zeros(NumElement*length(Line),1);
M=zeros(NumElement,numChannel);
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
%     ElementDensity(i)=calllib('libxrl','ElementDensity',Z(i));%
    density=ElementDensity(Z(i));%
%      A(i)=calllib('libxrl','AtomicWeight',Z(i));
    na(i)=density*NA/A(Z(i))*Z(i)*1e-24;
    j=1;
    while (j<=length(Line))
%         LineEnergy(i,j)=calllib('libxrl','LineEnergy',Z(i),-Line(j));
        new_energy=LineEnergy(Z(i),Line(j));%
        %         mu=na(i)*calllib('libxrl','CS_Total',Z(i),E0);
        %         Jump=calllib('libxrl','JumpFactor',Z(i),shell);
        %         w=calllib('libxrl','FluorYield',Z(i),shell);
        %         c=1; %%% weight fraction, 1: single element, ~1: compound
        %         intensity=I0*c*w*(1-1/Jump)*mu*T;
%                 CS_TotalBeam(i,1)=calllib('libxrl','CS_Total',Z(i),E0);
%         CS_FluoLine(i,j)=calllib('libxrl','CS_FluorLine',Z(i),-Line(j),E0);
        intensity=I0*CS_FluoLine(Z(i),Line(j));%
        BindingEnergy(length(Line)*(i-1)+j)=new_energy;
        [~,adj]=min(abs(repmat(new_energy,size(DetChannel))-DetChannel));
        PurePeak(adj)=PurePeak(adj)+intensity;
        
        j=j+1;
    end
    % fprintf('    %s           %d          %i          %5.2f               %5.2f\n', Element{Z(i)}, Z(i),-Line(j), new_energy, intensity);
    
    %%%%=============================== Start Gaussian Convolution
    M(i,:)=ifft(fft(PurePeak).*G);
    if(plotUnit)
        cmap=colormap(lines);
        cmap = cmap(1:NumElement,:);
        semilogy(DetChannel,abs(M(i,:)),'color',cmap(i,1:3),'LineStyle','-.','LineWidth',1.5);%
        Legend{i}=sprintf('%s',Element{Z(i)});
        legend(Legend)
        title('Gaussian spectrum for each element')
        hold on;
    end
    
    %%%%===================================================================
    %fprintf('=====================================================\n')
    i=i+1;
end
% for i=1:NumElement
%     for j=1:length(Line)
% 
%                 for t=1:NumElement
%                 CS_Total(i,j,t)=calllib('libxrl','CS_Total',Z(i),LineEnergy(t,j));
%                 end
%     end
% end
%  save(['./data/xRayLib',num2str(E0),'.mat'],'ElementDensity','LineEnergy','CS_FluoLine','CS_TotalBeam','CS_Total');
%  save('./data/AtomicWeight.mat','A');
 M=M.*1e-9;
if(TakeLog)
    M=abs(M);%1e-300;
end



