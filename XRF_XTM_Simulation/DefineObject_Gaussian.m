%%=======================================================================
%%% Input: Z: Atomic Numbers of Existing Elements in the object
%%%        E0: Beam Energy
%%%        m: Resolution of object
%%% Output: W: Concentration fraction of each Z in each pixel with
%%%            dimension m(1) by m(2) by cardinality(Z)
%%%         O: Object specified by W.*Z
%%%         MU: Attenuation matrix of O
%%=======================================================================
global x y dz m omega AbsorbScale MU_e Z 
global NumElement 
global slice onlyXRF

load PeriodicTable
if(synthetic)
    AbsorbScale=1; % scale on all attenuation coefficient
    ScaleM=1; % scale on unit spectra
else
    AbsorbScale=1;
    if(strcmp(sample,'Seed'))
        ScaleM=1;%4e1;
    else
        ScaleM=1;%5e3;
    end
    scale=0; % scale on fluorescence attenuation coefficient
end
%%%%%======================================
x=linspace(omega(1),omega(2),m(1)+1);
y=linspace(omega(3),omega(4),m(2)+1);
dz=[(omega(2)-omega(1))/m(2) (omega(4)-omega(3))/m(1)];
%%%========================== the grids of object
xc = getCellCenteredGrid(omega,[m(2) m(1)]);
xc = reshape(xc,prod(m),2);

%%%=========== assign weight matrix for each element in each pixel
if(synthetic)
    if(strcmp(sample,'Golosio') | strcmp(sample, 'checkboard'))
        Z = [6 8 14 20 26];% Golosio's Sample
    elseif(strcmp(sample,'Phantom'))
        Z = [19 31 26];
    elseif(strcmp(sample,'circle'))
        Z = [14];
    end
else
    if(strcmp(sample,'Rod'))
        pix_inner = find((xc(:,1)-center(2,1)).^2+(xc(:,2)-center(2,2)).^2 <= (1.5*r(2))^2);
        if(Si)
            Z=14;
        elseif(W_element)
            Z=74;
        else
            Z=[79 14 74];% Glass Rod
        end
    elseif(strcmp(sample,'Seed'))
        Z=[14 16 17 19 20 22 23 24 25 26 28 29 30 31 80 33 34 35 92 37 38 39 40];% Complete Seed
        Z=Z(slice-4);
    end
end
%---------------------------
NumElement=length(Z);
W=ones(N(level),N(level),NumElement);
if(level==1)
    if(synthetic)
        if(strcmp(sample,'Golosio'))
            CreateCircle; 
        elseif(strcmp(sample,'Phantom'))
            CreateElement; 
        elseif(strcmp(sample,'circle'))
            [X,Y]=meshgrid(1:m(1),1:m(2));
            center=[m(1)/2, m(2)/2];
            pix = (X-center(1)).^2+(Y-center(2)).^2 <= (m(1)/4)^2;
            % pix = (X-center(1)).^2+(Y-center(2)).^2 >= (m(1)/5)^2 & (X-center(1)).^2+(Y-center(2)).^2 <= (m(1)/4)^2;
            W(pix)=10;
        elseif(strcmp(sample,'checkboard'))
            W = kron(invhilb(N(1)/10)<0, ones(10,10));
            W=repmat(W,[1 1 NumElement]);
        end
        UnitSpectrumSherman_synthetic; %% Produce BindingEnergy M
    else
        for tsub=1:NumElement
            if(strcmp(sample,'Seed'))
                W(:,:,tsub)=abs(flipud(permute(iR_num(:,:,slice(tsub)),[2 1 3])));%tsub*2e-1;
            elseif(strcmp(sample,'Rod'))
                W(:,:,tsub)=abs(fliplr(rot90(permute(iR_num(:,:,tsub),[2 1 3]))));%tsub*2e-1;
            end
        end
        UnitSpectrumSherman_real; %% Produce BindingEnergy M
        clear iR_num iR
    end
    clear Line ElementDensity LineEnergy CS_FluoLine CS_Total CS_TotalBeam
end
%%%%%================== Attenuation Matrix at beam energy
MU_e=MU_e.*AbsorbScale; %% Discrete Scale
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,1,NumElement);
MU=sum(W.*repmat(MUe,[m(1),m(2),1]),3);
MU_XTM=MU;
%%%%% =================== Attenuation Matrix at flourescence energy (Corrected Attenuation)
MU_after=zeros(prod(m),NumElement);
for i=1:NumElement
    MU_after(:,i)=reshape(W,prod(m),NumElement)*MU_e(:,1,i+1);
end
%%%%% ====================================================================

%%=================== Picture the object based on atomic number and attenuation
if(PlotObject)
    figureObject(W,Z,m,NumElement,MU_e,0)
end

%%=======================================================================
% W(1,1,:)=[7.61/9 11/9 1.8/9 1.32/9 2.84/9 5/9 19/9];
% W(1,2,:)=[7.61/9 11/9 1.8/9 1.32/9 2.84/9 5/9 19/9];
%  W(1,3,:)=[7.61/9 11/9 1.8/9 1.32/9 2.84/9 5/9 19/9];
% W(2,1,:)= [7.61/9 11/9 1.8/9 1.32/9 2.84/9 5/9 19/9];
% W(2,2,:)=[7.61/9 11/9 1.8/9 1.32/9 2.84/9 5/9 19/9];
% W(2,3,:)=[7.61/9 11/9 1.8/9 1.32/9 2.84/9 5/9 19/9];
% W(3,1,:)=[7.61/9 11/9 1.8/9 1.32/9 2.84/9 5/9 19/9];
% W(3,2,:)=[7.61/9 11/9 1.8/9 1.32/9 2.84/9 5/9 19/9];
% W(3,3,:)=[7.61/9 11/9 1.8/9 1.32/9 2.84/9 5/9 19/9];
%%%========================================
% W(1,1,:)=[1 1 1 0.5 ];
% W(1,2,:)=[0 1 1 0.5];
% W(1,3,:)=[1 1 1 0];
% W(2,1,:)=[0 0.5 0.5 0.5];
% W(2,2,:)=[0.7 0.3 0 0.5];
% W(2,3,:)=[0.5 0.5 1 0];
% W(3,1,:)=[0.8 0 0 0.2];
% W(3,2,:)=[1 0 1 0.1];
% W(3,3,:)=[0 0.1 1 0];
%%%========================================
% W(1,1,:)=[1 0 0 0 ];
% W(1,2,:)=[0 1 0 0];
% W(1,3,:)=[1 0 0 0];
% W(1,4,:)=[1 0 0 0];
% W(2,1,:)=[0 0 0 1];
% W(2,2,:)=[0 0 1 0] ;
% W(2,3,:)=[0 1 0 0];
% W(2,4,:)=[0 1 0 0];

% W(3,1,:)=[0 0 1 0] ;
% W(3,2,:)= [1 0 0 0] ;
% W(3,3,:)=[0 1 0 0];
% W(3,4,:)=[0 1 0 0];
% W(4,1,:)=[0 0 1 0] ;
% W(4,2,:)= [1 0 0 0] ;
% W(4,3,:)=[0 1 0 0];
% W(4,4,:)=[0 1 0 0];
%%%========================== locate element attenuation coefficient
