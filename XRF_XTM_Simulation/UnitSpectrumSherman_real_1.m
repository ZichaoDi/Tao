
% Net intensity = I0*c_e*omega_e*(1-1/JumpFactor)*mu*T
% c_e: element concentration , for single element, c_e=1
% omega_e: fluorescence yield
% mu: linear attenuation coefficient at I0
% T: sample thickness, for single element, assume constant T=1 for all
% n_a: the number density of atoms
% mu=n_a*CrossSection
% K_shell:0 L1:1 L2:2 L3:3 M1:4 M2:5 M3:6 M4:7 M5:8
% KA_LINE:0 KB:1 LA:2 LB:3
global TakeLog M_decom M_raw Element synthetic
if(strcmp(sample,'Seed'))
    E0=12;
elseif(strcmp(sample,'Rod'))
    E0=12.1;
else
    E0=12.1;
end
E2I=1/3e8/6.62e-34*1e-15;
I0=E0*E2I;
NA=6.02e23;%Avogadro's number
% Z=1:95;
% NumElement=length(Z);
% numChannel_decom=length(Z);
load_xraylib=1;
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
% Line=[-219:-1];% all K,L,M lines
Line=[-3 -2 -1 -5 -6 -8 -13 -90 -34 -33 -102 -91 -98 -36 -35 -94 -89 -63 -95 -68 -207 -206];% partially K,L,M lines
% Line=[0:3];% K_alpha, K_beta, L_alpha, L_beta; %[-0 -1 -2 -3]; %% Transition Line, detailed defination see xraylib-lines.h
shell=[0  0  0  0  0  0   0   3   1   1    3   3   3   1   1   3   3   2   3   2  8    8];  %% Shell type
BindingEnergy=zeros(NumElement,length(Line));
M_decom=zeros(NumElement,numChannel_decom);
M_raw=zeros(NumElement,numChannel);
i=1;
T=1;
if(plotUnit)
    figure('name', 'Unit Spectrum');
end
%fprintf(1,'Element     AtomicNumber    Line     FluorescenceEnergy     Intensity\n')
%%==set up Gaussian Distribution
mu=0;
% FWHM = (DetChannel(2)-DetChannel(1))*20;
% sigma = FWHM/2.35;
truncInd_sub=cell(NumElement,1);
truncWidth=(DetChannel(2)-DetChannel(1))*5;% 0.02;
%%======================================================================
while(i<=NumElement)
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
            LineEnergy(Z(i),j)=calllib('libxrl','LineEnergy',Z(i),Line(j));
            % mu=na(i)*calllib('libxrl','CS_Total',Z(i),E0);
            % Jump=calllib('libxrl','JumpFactor',Z(i),shell(j));
            % w=calllib('libxrl','FluorYield',Z(i),shell(j));
            % c=1; %%% weight fraction, 1: single element, ~1: compound
            % CS_FluoLine(Z(i),j)=c*w*(1-1/Jump)*mu*T;
            CS_FluoLine(Z(i),j)=calllib('libxrl','CS_FluorLine_Kissel',Z(i),Line(j),E0);
        end
        new_energy=LineEnergy(Z(i),j);%
        intensity=CS_FluoLine(Z(i),j);%
        BindingEnergy(i,j)=new_energy;
        [~,adj]=sort(abs(repmat(new_energy,size(DetChannel))-DetChannel*1e3));
        PurePeak=0.*DetChannel';
        PurePeak(adj(1))=intensity;    
        M_decom(i,i)=M_decom(i,i)+intensity;
        res=abs(DetChannel(adj(2))-DetChannel(adj(1))); % KeV
        sigma=sqrt((0.1/2.3584)^2+3.58*1e-3*0.114*BindingEnergy(i,j));
        G=1/(sigma*sqrt(2*pi))*exp(-(DetChannel'-BindingEnergy(i,j)).^2./(2*sigma^2));
        M_raw(i,:)=M_raw(i,:)+ifft(fft(PurePeak).*fft(G));
        truncInd_sub{i}=unique([truncInd_sub{i}; find(DetChannel> BindingEnergy(i,j)-truncWidth & DetChannel < BindingEnergy(i,j)+truncWidth & DetChannel <= E0)]);
        j=j+1;
    end
    truncInd_sub{i}=sort(truncInd_sub{i});
    
    %%%%=============================== Start Gaussian Convolution
    if(plotUnit)
        cmap=colormap(lines);
        cmap = cmap(1:NumElement,:);
        plot(DetChannel,abs(M_raw(i,:)),'color',cmap(i,1:3),'LineStyle','-.','LineWidth',1.5);
        Legend{i}=sprintf('%s',Element{Z(i)});
        legend(Legend)
        title('Gaussian spectrum for each element')
        hold on;
    end
    
    %fprintf('=====================================================\n')
    i=i+1;
end
truncInd=[];
for i=1:NumElement
truncInd=[truncInd;truncInd_sub{i}];
end
truncInd=sort(unique(truncInd));
%%%%===================================================================
if(load_xraylib)
    for i=1:NumElement
     for j=1:length(Line)
         for t=1:NumElement
             CS_Total(Z(i),j,Z(t))=calllib('libxrl','CS_Total',Z(i),LineEnergy(Z(t),j));
         end
     end
    end
    % save(['./data/xRayLib',num2str(E0),'.mat'],'ElementDensity','LineEnergy','CS_FluoLine','CS_TotalBeam','CS_Total');
end
%%=======================================================================
MU_e=zeros(NumElement,length(Line),1+NumElement);
for i=1: NumElement
    if(load_xraylib)
        MU_e(i,1,1)=calllib('libxrl','CS_Total',Z(i),E0);
    else
        MU_e(i,1,1)=CS_TotalBeam(Z(i),1);% units: cm^2/g
    end
    for j=1:NumElement
        for j_line = 1:length(Line)
            if(load_xraylib)
                MU_e(i,j_line,j+1)=calllib('libxrl','CS_Total',Z(i),BindingEnergy(j,j_line));% only consider K_alpha line
            else
                MU_e(i,j_line,j+1)=CS_Total(Z(i),j_line,Z(j));%
            end
        end
    end
end
