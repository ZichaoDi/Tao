% function [O,MU,W]=DefineObject(Z,E0)
%%=======================================================================
%%% Input: Z: Atomic Numbers of Existing Elements in the object
%%%        E0: Beam Energy
%%%        m: Resolution of object
%%% Output: W: Concentration fraction of each Z in each pixel with
%%%            dimension m(1) by m(2) by cardinality(Z)
%%%         O: Object specified by W.*Z
%%%         MU: Attenuation matrix of O
%%=======================================================================
global x y m omega dz AbsorbScale min_MU max_MU
global XTMscale NumLines NumElement 
global whichElement slice onlyXRF


AbsorbScale=1e-0;
ScaleM=1e-5;
%%%%%======================================
DisY=m(1)+1;
DisX=m(2)+1;
x=linspace(omega(1),omega(2),DisX);
y=linspace(omega(3),omega(4),DisY);
dz=[(omega(2)-omega(1))/m(2) (omega(4)-omega(3))/m(1)];
[X,Y] = meshgrid((x(1:end-1)+x(2:end))./2,(y(1:end-1)+y(2:end))./2);
center=[0 0];
%%%========================== the grids of object
xc = getNodalGrid(omega,[m(2) m(1)]);
%%%========================== assign weight matrix for each element in each pixel
if(onlyXRF)
    if (xrf_roi)
        NumElement=1;
        data=Phantom3_Young(:,:,:,whichElement+1);
        W=data(1:current_n,1:current_n,slice);
    else
        NumElement=4;
        W=zeros(current_n,current_n,NumElement);
        for i=1:NumElement
            data=Phantom3_Young(:,:,:,i+1);
            W(:,:,i)=data(1:current_n,1:current_n,slice);
        end
    end
else
    CreateElement; %% shepp-logan phantom
%     CreateCircle; %% sample with circles 
    %------------------------- a sample to test the different impact from heavy and light elements
    % SvenSample;
    %----------------------------
    % load W_sample10
    % W=W_sample10;
    %---------------------------
%     NumElement=2;
%     W=zeros(m(1),m(2),NumElement);
%     for tsub=1:1;%NumElement
%         W(:,:,tsub)=tsub*2e-1;
%     end
end
% NumElement=2;
%%%%%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Z=[29 30 46 49 57 74 79];%% 42 29 26 ];%20 49 57 46];% reference sample: Pb La Pd Mo Cu Fe Ca
% ComChoices=nchoosek(1:6,3);
% Z=Z(ComChoices(1,:));
Z=Z(1:NumElement);
UnitSpectrumSherman_Gaussian; %% Produce BindingEnergy M
%%=======================================================================
NumLines=NumElement;
MU_e=zeros(NumElement,1,1+NumLines);
for i=1: NumLines
    MU_e(i,1,1)=na(i)*CS_TotalBeam(Z(i),1);%calllib('libxrl','CS_Total',Z(i),E0);
    for j=1:NumElement
        MU_e(i,1,j+1)=na(i)*CS_Total(Z(i),1,Z(j));%calllib('libxrl','CS_Total',Z(i),BindingEnergy(j));   only consider K_alpha line
    end
end
%%%%%================== Attenuation Matrix at beam energy
MU_e=MU_e.*AbsorbScale; %% Discrete Scale
min_MU=min(MU_e(:,1,1))-1;%=-1./log(MU_e(:,1,1));
max_MU=max(MU_e(:,1,1));
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(MU_e(:,1,1),1,1,NumElement);
MU=sum(W.*repmat(MUe,[m(1),m(2),1]),3);
XTMscale=1e0;
MU_XTM=MU.*XTMscale;
%%%%% =================== Attenuation Matrix at flourescence energy (Corrected Attenuation)
MU_after=cell(NumElement,1);
for i=1:NumElement
    MU_after{i}=sum(W.*repmat(reshape(MU_e(:,1,i+1),1,1,NumElement),[m(1),m(2),1]),3);
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