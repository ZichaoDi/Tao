% Net intensity = I0*c_e*omega_e*(1-1/JumpFactor)*mu*T
% c_e: element concentration , for single element, c_e=1
% omega_e: fluorescence yield
% mu: linear attenuation coefficient at I0
% T: sample thickness, for single element, assume constant T=1 for all
% n_a: the number density of atoms
% mu=n_a*CrossSection
% K_shell:0 L1:1 L2:2 L3:3 M1:4 M2:5 M3:6 M4:7 M5:8
% KA_LINE:0 KB:1 LA:2 LB:3
global TakeLog I0 M_decom M_raw Element synthetic
numChannel_decom=NumElement;
DetChannel_decom=[1:NumElement]';
if(synthetic)
E0=12.1;
E2I=1/3e8/6.62e-34*1e-15;
I0=E0*E2I;
else
E0=12.1;
I0=12.1;
end
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
% Line=0:3;%[-0 -1 -2 -3]; %% Transition Line, detailed defination see xraylib-lines.h
% Line=-[0:207];
Line=[-3 -2 -1 -5 -6 -8 -13 -90 -34 -33 -102 -91 -98 -36 -35 -94 -89 -63 -95 -68 -207 -206];% partially K,L,M lines

shell=0;  %% Shell type
BindingEnergy=zeros(NumElement,length(Line));
M_decom=zeros(NumElement,numChannel_decom);
i=1;
T=1;
%fprintf(1,'Element     AtomicNumber    Line     FluorescenceEnergy     Intensity\n')
%%==set up Gaussian Distribution
mu=0;
%%======================================================================
while(i<=NumElement)
PurePeak=0.*[1:NumElement]';
if(load_xraylib)
ElementDensity(Z(i))=calllib('libxrl','ElementDensity',Z(i));%
A(Z(i))=calllib('libxrl','AtomicWeight',Z(i));
CS_TotalBeam(Z(i),1)=calllib('libxrl','CS_Total',Z(i),E0);
end
density=ElementDensity(Z(i));%
j=1;
while (j<=length(Line))
if(load_xraylib)
LineEnergy(Z(i),j)=calllib('libxrl','LineEnergy',Z(i),Line(j));
CS_FluoLine(Z(i),j)=calllib('libxrl','CS_FluorLine_Kissel',Z(i),Line(j),E0);
end
new_energy=LineEnergy(Z(i),j);%
intensity=CS_FluoLine(Z(i),j);%
BindingEnergy(i,j)=new_energy;
[~,adj]=min(abs(repmat(new_energy,size(DetChannel_decom))-DetChannel_decom));
PurePeak(adj)=PurePeak(adj)+intensity;    
M_decom(i,i)=M_decom(i,i)+sum(PurePeak);
j=j+1;
end
%fprintf('=====================================================\n')
i=i+1;
end
%%%%===================================================================
if(load_xraylib)
    for i=1:NumElement
        for j=1:length(Line)
            for t=1:NumElement
                 CS_Total(Z(i),j,Z(t))=calllib('libxrl','CS_Total',Z(i),LineEnergy(Z(t),j));
            end
        end
    end
end
%%=======================================================================
NumLines=NumElement;
MU_e=zeros(NumElement,1,1+NumLines);
for i=1: NumLines
    if(load_xraylib)
        MU_e(i,1,1)=calllib('libxrl','CS_Total',Z(i),E0);
    else
        MU_e(i,1,1)=CS_TotalBeam(Z(i),1);% units: cm^2/g
    end
    for j=1:NumElement
        if(load_xraylib)
            MU_e(i,1,j+1)=calllib('libxrl','CS_Total',Z(i),BindingEnergy(j,1));% only consider K_alpha line
        else
            MU_e(i,1,j+1)=CS_Total(Z(i),1,Z(j));%
        end
    end
end
