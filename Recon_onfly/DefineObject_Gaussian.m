%%=======================================================================
%%% Input: Z: Atomic Numbers of Existing Elements in the object
%%%        E0: Beam Energy
%%%        m: Resolution of object
%%% Output: W: Concentration fraction of each Z in each pixel with
%%%            dimension m(1) by m(2) by cardinality(Z)
%%%         O: Object specified by W.*Z
%%%         MU: Attenuation matrix of O
%%=======================================================================
global xc x_coor y_coor dz m omega 
global MU_e Z 
global NumElement Element 
global slice onlyXRF

load PeriodicTable
%%%%%======================================
x_coor=linspace(omega(1),omega(2),m(1)+1);
y_coor=linspace(omega(3),omega(4),m(2)+1);
%%%========================== the grids of object
xc = getCellCenteredGrid(omega,[m(2) m(1)]);
xc = reshape(xc,prod(m),2);

%%%=========== assign weight matrix for each element in each pixel
if(synthetic)
    if(strcmp(sample,'Golosio') | strcmp(sample, 'checkboard'))
        Z = [6 8 14 20 26];% Golosio's Sample
    elseif(strcmp(sample,'Phantom'))
        Z = [19];% 31 26];
    elseif(strcmp(sample,'circle'))
        Z = [14];
    elseif(strcmp(sample, 'fakeRod'))
        NumElement=3;
        Z=[79 5 8 11 13 14 19 74];% Glass Rod
        if(NumElement==3)
            Z=Z([1 6 8]);
        end
    end
else
    if(strcmp(sample,'Rod'))
        Z=[79 14 74];% Glass Rod
    elseif(strcmp(sample,'Seed'))
        Z=[14 16 17 19 20 22 23 24 25 26 28 29 30 31 80 33 34 35 92 37 38 39 40];% Complete Seed
        Z=Z(slice-4);
    end
end
%---------------------------
NumElement=length(Z);
    W=ones(N(1),N(1),NumElement);
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
        elseif(strcmp(sample,'fakeRod'))
            CreateRod;
        end
    else
        for tsub=1:NumElement
            if(strcmp(sample,'Seed'))
                W(:,:,tsub)=abs(fliplr(permute(iR_num(:,:,tsub),[2 1 3])));%tsub*2e-1;
            elseif(strcmp(sample,'Rod'))
                W(:,:,tsub)=abs(fliplr(rot90(permute(iR_num(:,:,tsub),[2 1 3]))));%tsub*2e-1;
            end
        end
    end
    if(onlyXRF)
        UnitSpectrumSherman_real_1; %% Produce BindingEnergy M
    else
        UnitSpectrumSherman_real; %% Produce BindingEnergy M
    end
    clear Line ElementDensity LineEnergy CS_FluoLine CS_Total CS_TotalBeam
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,1,NumElement);
MU=sum(W.*repmat(MUe,[m(1),m(2),1]),3);
MU_XTM=MU;
%%%%% =================== Attenuation Matrix at flourescence energy (Corrected Attenuation)
MU_after=zeros(prod(m),NumElement);
for i=1:NumElement
    temp=sum(reshape(W,prod(m),NumElement)*MU_e(:,:,i+1),2);
    temp=flipud(reshape(temp,m(1),m(2))');
    MU_after(:,i)=temp(:);
end
clear temp;
%%%%% ====================================================================

%%=================== Picture the object based on atomic number and attenuation
if(PlotObject)
    figureObject(W,Z,m,NumElement,MU_e,0)
end

