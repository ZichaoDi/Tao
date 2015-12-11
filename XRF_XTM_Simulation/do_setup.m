%---------------------------------------------------------------
% Set up a problem for optimization [multigrid or regular].
%---------------------------------------------------------------
global N numThetan NF
global bounds LogScale Joint
global grad_type err0 WS W0
global onlyXRF NoSelfAbsorption
global W_level xrf_level xtm_level L_level GI_level SI_level SigmaR SigmaT m_level nTau_level
global NumElement xinitial current_n
%--------------------------------------------------
% Select technique for gradient calculation.

grad_type = 'adj';  % 'adj' = adjoint/exact
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
end
N=20;% [33 17 9 5 3];% 17 9];%[129 65  9 5];% 
if(synthetic)
    numThetan=5; % number of scanning angles/projections
    numChannel=500;% number of energy channels on xrf detector
else
    if(strcmp(sample,'Seed'))
        load ./data/ApsDataExtract/DogaSeeds/DownSampleSeeds221_elements.mat
        data = permute(data_H,[2,3,1]);
        slice = 5:27;  
        data_h=[];
        ang_rate=30;
        tau_rate=10;
        for ele=1:size(data,1)
            data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
        end
        data_xrt=data_h(3,:,:); %transimission data
        data_xrf=data_h(slice,:,:);
        DecomposedElement=1;
        [x_ir,y_ir]=meshgrid(1:size(iR,1));
        [x_num,y_num]=meshgrid(linspace(1,size(iR,1),N(1)));
        iR_num=zeros(N(1),N(1),size(iR,3));
        for ele=1:size(iR,3)
            iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,ele),x_num,y_num);
        end
        clear x_ir y_ir x_num y_num
    elseif(strcmp(sample,'Rod'))
        % load ~/Documents/Research/APSdata/GlassRod/2dSlice/Slice1
        % load ~/Documents/Research/APSdata/GlassRod/2dSlice/rec_XTM_59
        % iR=imread('~/Documents/Research/APSdata/GlassRod/Au_L_3/Au_L_0000.tif');
        slice = [25 10 32 30 3 34]; %GlassRod% 
    end
    numThetan=size(data_h,2);
end
%%-----------------------------------------------
%%-----------------------------------------------
NoSelfAbsorption=0;
nh=N(1); 
n_level=1;
level=[1:n_level];
N=nh./(2.^(level-1))+(2.^(level-1)-1)./2.^(level-1);
NF = [0*N; 0*N; 0*N];
nm=length(N);
angleScale=1; %1: half angle; 2: full angle
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
bounds = 1;  % no bound constraints
Joint=1; % 0: XRF; -1: XTM; 1: Joint inversion
ReconAttenu = 1; % 0: Recover W; 1: Recover miu
Alternate=1;
LogScale=1; %% determine if the XTM is solved taking log first or not
Weighted=0; %% 1 if use weighted least-square form
linear_S=0;
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
    if(level==1)
        XRF_XTM_Tensor;
        WS=W;%MU;% 
    else
        W=downdate(W(:),1);
        W=reshape(W,N(level),N(level),NumElement);
    end
    % XTM_Tensor;
    W_level{level}=W;
    xtm_level{level}=DisR;
    
    L_level{level}=L;
    GI_level{level}=GlobalInd;
    nTau_level(level)= nTau;
    m_level(level,:)=m;
    SigmaT{level}=SigMa_XTM;
%-----------------------------------------
    xrf_level{level}=XRF;
    SI_level{level}=SelfInd;
    SigmaR{level}=SigMa_XRF;
end

nTol=N(1)^2*NumElement;
%---------------------------------------------
% Specify initial guess for optimization.
% xinitial=x0;
