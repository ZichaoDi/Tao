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
if(synthetic)
    if(strcmp(sample,'Golosio') | strcmp(sample, 'checkboard'))
        Z = 19;%[6 8 14 20 26];% Golosio's Sample
    elseif(strcmp(sample,'Phantom')|strcmp(sample,'mri'))
        Z = [19];
    elseif(strcmp(sample,'circle'))
        Z = [14];
    elseif(strcmp(sample, 'fakeRod'))
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
    elseif(strcmp(sample,'Paunesku'))
        Z=[15 16 20 22 26 30];
        Z=Z(ele_ind);
    elseif(strcmp(sample,'Filter'))
        Z=[14 15 16 20 30 58];
    elseif(strcmp(sample,'miller'))
        Z=[20 30 68];
    elseif(strcmp(sample,'Run02'))
        Z=20;
    elseif(strcmp(sample,'20173'))
        Z=[15 16 17 19 20 25 26 29 30];
        Z=Z(slice_tot);
    elseif(strcmp(sample,'Zn_modified_6'))
        Z=30;
    end
end
%---------------------------
NumElement=length(Z);
iR_num=zeros(N,N,NumElement);
    W=zeros(N(1),N(1),NumElement);
    if(synthetic)
        if(strcmp(sample,'Golosio'))
            CreateCircle; 
        elseif(strcmp(sample,'Phantom')|strcmp(sample,'mri'))
            CreateElement; 
        elseif(strcmp(sample,'circle'))
            [X,Y]=meshgrid(1:m(1),1:m(2));
            center=[m(1)/3, m(2)/3];
            r=5;
            pix = (X-center(1)).^2+(Y-center(2)).^2 <= r^2;%& (X-center(1)).^2+(Y-center(2)).^2>=(r-5).^2; %% circle
            W(pix)=1;
        elseif(strcmp(sample,'checkboard'))
            W = kron(invhilb(round(N(1)/10))<0, ones(10,10));
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
            elseif(strcmp(sample,'Filter'))
                W(:,:,tsub)=abs((rot90(iR_num(:,:,tsub),3)));
    elseif(strcmp(sample,'Paunesku'))
                W(:,:,tsub)=abs(fliplr(rot90(permute(iR_num(:,:,tsub),[2 1 3]))));%tsub*2e-1;
            end
        end
        clear iR
    end
    UnitSpectrumSherman_real; %% Produce BindingEnergy M
    clear Line ElementDensity LineEnergy CS_FluoLine CS_Total CS_TotalBeam
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,1,NumElement);
MU_XTM=sum(W.*repmat(MUe,[m(1),m(2),1]),3);
%%%%% ====================================================================

%%=================== Picture the object based on atomic number and attenuation

