%---------------------------------------------------------------
% Set up a problem for optimization [multigrid or regular].
%---------------------------------------------------------------
global N numThetan NF
global bounds LogScale Joint
global grad_type err0 WS W0
global onlyXRF NoSelfAbsorption
global W_level xrf_level xtm_level L_level GI_level SI_level SigmaR SigmaT m_level nTau_level
global NumElement xinitial current_n Tik penalty
%--------------------------------------------------
% Select technique for gradient calculation.

grad_type = 'full-linear';  % 'adj' = adjoint/exact
% 'fdr' = finite-difference [real]
% 'fdi' = finite-difference [complex]

%--------------------------------------------
% Initialize arrays for discretizations
Tomo_startup;
%%===============Load Sample=====================
synthetic=0;
if(synthetic)
    sample='Phantom';
else
    sample='Seed';
    % sample='Rod';
end
N=70;% [33 17 9 5 3];% 17 9];%[129 65  9 5];%
angleScale=2; %1: half angle; 2: full angle
if(synthetic)
    numThetan=15; % number of scanning angles/projections
    numChannel=5;% number of energy channels on xrf detector
    DecomposedElement=0;
else
    if(strcmp(sample,'Seed'))
        load ./data/ApsDataExtract/DogaSeeds/DownSampleSeeds441_elements.mat
        load RawSpectra_seed.mat
        if(angleScale==1)
            cutInd=floor(size(data_H,3)/2);
            data_H=data_H(:,:,1:cutInd);
            spectra=spectra(:,:,1:cutInd);
        end
        data = permute(data_H,[2,3,1]);
        slice =[13 14 17];%5:27;%      
        data_h=[];
        ang_rate=10;
        tau_rate=6;
        for ele=1:size(data,1)
            data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
        end
        if(ndims(data_h)==2)
            data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
        end
        DecomposedElement=1;
        truncChannel=(DecomposedElement==0)*1;
        data_xrt=data_h(3,:,:); %transimission data
        if(DecomposedElement)
            data_xrf=data_h(slice,:,:);
        else
            data_xrf=permute(spectra(1:tau_rate:end,:,1:ang_rate:end),[2 3 1]);
        end
        [x_ir,y_ir]=meshgrid(1:size(iR,1));
        [x_num,y_num]=meshgrid(linspace(1,size(iR,1),N(1)));
        iR_num=zeros(N(1),N(1),size(iR,3));
        for ele=1:size(iR,3)
            iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,ele),x_num,y_num);
        end
    elseif(strcmp(sample,'Rod'))
        load ~/Documents/Research/APSdata/GlassRod/2dSlice/Slice30
        load ~/Documents/Research/APSdata/GlassRod/2dSlice/tomoRec_3
        slice_tot = [3 4 25 30]; %GlassRod%
        slice = [4 25 30];
        data_h=[];
        ang_rate=10;
        tau_rate=45;
        for ele=1:size(data,1)
            data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
        end
        if(ndims(data_h)==2)
            data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
        end
        DecomposedElement=1;
        truncChannel=(DecomposedElement==0)*0;
        data_xrt=data_h(28,:,:); %transimission data: V_dpc_cfg
        data_xrf=[];
        if(DecomposedElement)
            data_xrf(1,:,:)=data_h(slice_tot(1),:,:)+data_h(slice_tot(2),:,:);
            data_xrf(2,:,:)=data_h(slice_tot(3),:,:);
            data_xrf(3,:,:)=data_h(slice_tot(4),:,:);
        else
            load(['data/ApsDataExtract/2xfm1211_14/2xfm_',num2str(172+1),'.mat'],'XRF','DetChannel');
            data_xrf=permute(XRF(:,30,:),[3 2 1]);%permute(spectra(1:tau_rate:end,:,1:ang_rate:end),[2 3 1]);
            clear XRF
        end
        [x_ir,y_ir]=meshgrid(1:size(iR,1));
        [x_num,y_num]=meshgrid(linspace(1,size(iR,1),N(1)));
        iR_num=zeros(N(1),N(1),size(iR,3));
        for ele=1:size(iR,3)
            iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,ele),x_num,y_num);
        end
    end
    clear x_ir y_ir x_num y_num iR data spectra 
    numThetan=size(data_h,2);
end
%%-----------------------------------------------
nh=N(1);
n_level=1;
level=[1:n_level];
N=nh./(2.^(level-1))+(2.^(level-1)-1)./2.^(level-1);
NF = [0*N; 0*N; 0*N];
nm=length(N);
W_level=cell(nm,1);
xrf_level=cell(nm,1);
xtm_level=cell(nm,1);
L_level=cell(nm,1);
GI_level=cell(nm,1);
SI_level=cell(nm,1);
SigmaR=cell(nm,1);
SigmaT=cell(nm,1);
m_level=zeros(nm,2);
nTau_level= zeros(nm,1);
%%=============================
NoSelfAbsorption=0; % 0: include self-absorption in the XRF inversion
bounds = 1;  % no bound constraints
Joint=1; % 0: XRF; -1: XTM; 1: Joint inversion
ReconAttenu = 1*(Joint==-1); % 0: Recover W; 1: Recover miu
Alternate=1*(Joint~=-1);
linear_S=0*Alternate;
LogScale=1; %% determine if the XTM is solved taking log first or not
Weighted=0; %% 1 if use weighted least-square form
%--------------------------------------------------------------------
% Compute the dependent-variable arrays
PlotObject=0;
plotSpec = 0; % Do you want to see the spectra? If so plotSpec = 1
plotTravel=0; % If plot the intersection of beam with object
plotUnit=0;
plotElement=0;
plotResult=1;
%%------------------------------ Downsample data from finest level
% for level=1:nm
%     current_n=N(level);
%
%     if(level==1)
%         XTM_Tensor;
%         WS=MU;%W;
%         xtm_level{level}=DisR;
%         L_level{level}=L;
%         GI_level{level}=GlobalInd;
%         nTau_level(level)= nTau;
%     else
%         xtm_level{level}=downdate(xtm_level{level-1},3);
%         nTau_level(level)=size(xtm_level{level},1)-1;
%         XTM_Tensor;
%         L_level{level}=L;
%         GI_level{level}=GlobalInd;
%     end
%     W_level{level}=W;
%
%     m_level(level,:)=m;
%     SigmaT{level}=SigMa_XTM;
%     %-----------------------------------------
% %     xrf_level{level}=XRF;
% %     SI_level{level}=SelfInd;
% %     SigmaR{level}=SigMa_XRF;
% end
%%------------------------------ Use same finest data for each level
for level=1:nm
    current_n=N(level);
    if(Joint==-1)
        if(level==1)
            % [L,dd,~,~]=build_weight_matrix(MU,thetan,1,'area');
            XTM_Tensor;
            Lap{level}=delsq(numgrid('S',N(level)+2));
            L_level{level}=L;
        else
            IhH{level} = downdate_L(N, level);
            Lap{level}=IhH{level}*Lap{level-1}*IhH{level}';
            L_level{level}=L_level{level-1}*IhH{level}';
        end
        if(ReconAttenu)
            NumElement=1;
            if(level==1)
                W_level{level}=MU;
                WS= W_level{1};
            else
                W_level{level}=N(level)^2;
            end
        else
            W_level{level}=W;
            if(level==1)
                WS= W;
            end
        end
    else
        XRF_XTM_Tensor;
        W_level{level}=W;
        L_level{level}=L;
        if(level==1)
            WS=W;
        end
    end
    xtm_level{level}=DisR;
    GI_level{level}=GlobalInd;
    nTau_level(level)= nTau;
    m_level(level,:)=[N(level),N(level)];
    SigmaT{level}=SigMa_XTM;
    %-----------------------------------------
    if(Joint~=-1)
        xrf_level{level}=XRF;
        SI_level{level}=SelfInd;
        SigmaR{level}=SigMa_XRF;
    end
end
clear L W DisR GlobalInd SigMa_XTM SelfInd XRF SigMa_XRF data_h;
%%%================= First-order derivative regularization
penalty=0;
if(penalty)
    Tik=spalloc(length(W(:))-1,length(W(:)),2*(length(W(:))-1)); 
    for i=1:length(W(:))-1
        Tik(i,i)=-1; Tik(i,i+1)=1;
        % Tik(i,i)=-1;Tik(i,i+1)=2; Tik(i,i+2)=-1;
    end
end
nTol=N(1)^2*NumElement;
save([sample,num2str(N(1)),'_',num2str(numThetan),'_',num2str(nTau+1),'_',num2str(numChannel),'.mat'],'-v7.3')
