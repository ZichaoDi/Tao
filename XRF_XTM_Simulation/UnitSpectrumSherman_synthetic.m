
% Net intensity = I0*c_e*omega_e*(1-1/JumpFactor)*mu*T
% c_e: element concentration , for single element, c_e=1
% omega_e: fluorescence yield
% mu: linear attenuation coefficient at I0
% T: sample thickness, for single element, assume constant T=1 for all
% n_a: the number density of atoms
% mu=n_a*CrossSection
% K_shell:0 L1:1 L2:2 L3:3 M1:4 M2:5 M3:6 M4:7 M5:8
% KA_LINE:0 KB:1 LA:2 LB:3
global TakeLog I0 M Element synthetic

E0=12.1;
E2I=1/3e8/6.62e-34*1e-15;
I0=E0*E2I;
NumElement=length(Z);
NA=6.02e23;%Avogadro's number
load_xraylib=0;
if(load_xraylib)
    if(ismac)
         loadlibrary('/opt/local/lib/libxrl.dylib','/opt/local/include/xraylib/xraylib.h');
    else
         loadlibrary('/homes/wendydi/local/xraylib/lib/libxrl.so','/homes/wendydi/local/xraylib/include/xraylib/xraylib.h');
    end
else
    load(['xRayLib',num2str(E0),'.mat'])
end
load AtomicWeight
Line=[-3 -2 -6 -90 -89 -63 -95 -68 -207];%0:3;%[-0 -1 -2 -3]; %% Transition Line, detailed defination see xraylib-lines.h
shell=0;  %% Shell type
BindingEnergy=zeros(NumElement,length(Line));
M=zeros(NumElement,numChannel);
i=1;
T=1;
if(plotUnit)
    figure('name', 'Unit Spectrum');
end
%fprintf(1,'Element     AtomicNumber    Line     FluorescenceEnergy     Intensity\n')
%%==set up Gaussian Distribution
mu=0;
sigma=(DetChannel(2)-DetChannel(1))/2.35*20;
%%======================================================================
while(i<=NumElement)
    PurePeak=0.*DetChannel;
    if(load_xraylib)
        ElementDensity(Z(i))=calllib('libxrl','ElementDensity',Z(i));%
        A(Z(i))=calllib('libxrl','AtomicWeight',Z(i));
        CS_TotalBeam(Z(i),1)=calllib('libxrl','CS_Total',Z(i),E0);
    end
    density=ElementDensity(Z(i));%
    na(i)=density*NA/A(Z(i))*Z(i)*1e-24;
    j=1;
    while (j<=length(Line))
        if(load_xraylib)
            LineEnergy(Z(i),Line(j)+1)=calllib('libxrl','LineEnergy',Z(i),Line(j));
            CS_FluoLine(Z(i),Line(j)+1)=calllib('libxrl','CS_FluorLine',Z(i),Line(j),E0);
        end
        new_energy=LineEnergy(Z(i),Line(j)+1);%
        %         mu=na(i)*calllib('libxrl','CS_Total',Z(i),E0);
        %         Jump=calllib('libxrl','JumpFactor',Z(i),shell);
        %         w=calllib('libxrl','FluorYield',Z(i),shell);
        %         c=1; %%% weight fraction, 1: single element, ~1: compound
        %         intensity=I0*c*w*(1-1/Jump)*mu*T;
        intensity=I0*CS_FluoLine(Z(i),Line(j)+1);%
        BindingEnergy(i,j)=new_energy;
        [~,adj]=min(abs(repmat(new_energy,size(DetChannel))-DetChannel));
        PurePeak(adj)=PurePeak(adj)+intensity;    
        M(i,:)=M(i,:)+intensity/(sigma*sqrt(2*pi))*exp(-(DetChannel'-BindingEnergy(i,j)).^2./(2*sigma^2));;
        j=j+1;
    end
    % fprintf('    %s           %d          %i          %5.2f               %5.2f\n', Element{Z(i)}, Z(i),-Line(j), new_energy, intensity);
    
    %%%%=============================== Start Gaussian Convolution
    if(plotUnit)
        cmap=colormap(lines);
        cmap = cmap(1:NumElement,:);
        % semilogy(DetChannel,abs(M(i,:)),'color',cmap(i,1:3),'LineStyle','-.','LineWidth',1.5);
        plot(DetChannel,abs(M(i,:)),'color',cmap(i,1:3),'LineStyle','-.','LineWidth',1.5);
        Legend{i}=sprintf('%s',Element{Z(i)});
        legend(Legend)
        title('Gaussian spectrum for each element')
        hold on;
    end
    
    %fprintf('=====================================================\n')
    i=i+1;
end
%%%%===================================================================
if(load_xraylib)
    for i=1:NumElement
     for j=1:length(Line)
         for t=1:NumElement
             CS_Total(Z(i),Line(j)+1,Z(t))=calllib('libxrl','CS_Total',Z(i),LineEnergy(Z(t),Line(j)+1));
         end
     end
    end
    save(['./data/xRayLib',num2str(E0),'.mat'],'ElementDensity','LineEnergy','CS_FluoLine','CS_TotalBeam','CS_Total');
end
%%=======================================================================
NumLines=NumElement;
MU_e=zeros(NumElement,1,1+NumLines);
for i=1: NumLines
    if(load_xraylib)
        MU_e(i,1,1)=calllib('libxrl','CS_Total',Z(i),E0);% na(i)*
    else
        MU_e(i,1,1)=CS_TotalBeam(Z(i),1);%
    end
    for j=1:NumElement
        if(load_xraylib)
            MU_e(i,1,j+1)=calllib('libxrl','CS_Total',Z(i),BindingEnergy(j,1));% na(i)*   only consider K_alpha line
        else
            MU_e(i,1,j+1)=CS_Total(Z(i),1,Z(j));%
        end
    end
end
M=M.*ScaleM;
if(TakeLog)
    M=abs(M);%1e-300;
end

