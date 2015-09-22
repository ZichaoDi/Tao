%---------------------------------------------------------------
% Set up a problem for optimization [multigrid or regular].
%---------------------------------------------------------------
% W,xrfData,xtmData,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0

global N numThetan NF
global bounds LogScale Joint ReconAttenu
global grad_type err0 WS Beta W0 
global NoSelfAbsorption synthetic
global W_level xrf_level xtm_level L_level GI_level SI_level SigmaR SigmaT m_level nTau_level
global DecomposedElement NumElement xinitial current_n
%--------------------------------------------------
% Select technique for gradient calculation.

grad_type = 'adj';  % 'adj' = adjoint/exact
% 'fdr' = finite-difference [real]
% 'fdi' = finite-difference [complex]
solver = 'TN'; % 'GN'= Gauss-Newton; 
%--------------------------------------------
% Initialize arrays for discretizations
Tomo_startup;
%%-----------------------------------------------
load ./data/ApsDataExtract/DogaSeeds/DownSampleSeeds221_elements.mat
synthetic=1;
if(synthetic)
    sample='golosio';
    data_H=ones(221,27,725);
else
    sample='seed';
end
data = permute(data_H,[2,3,1]);
%%-----------------------------------------------
% load ~/Documents/Research/APSdata/GlassRod/2dSlice/Slice1
% load ~/Documents/Research/APSdata/GlassRod/2dSlice/rec_XTM_59
% iR=imread('~/Documents/Research/APSdata/GlassRod/Au_L_3/Au_L_0000.tif');
%%-----------------------------------------------
data_h=[];
ang_rate=10;
tau_rate=3;
for ele=1:size(data,1)
data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
end
slice = 5:27;% Seeds [25 10 32 30 3 34]; %GlassRod%  
data_xrt=data_h(3,:,:); %transimission data
data_xrf=data_h(slice,:,:);
DecomposedElement=1;
% data=h5read('~/Documents/MATLAB/APSdata/xfm_Doga/xfm_data_elements.h5','/exchange/data');
% data=squeeze(sum(data,2));
%%-----------------------------------------------
NoSelfAbsorption=1;
N= [65 33];% 9];% 33 17 9];% 
[x_ir,y_ir]=meshgrid(1:size(iR,1));
[x_num,y_num]=meshgrid(linspace(1,size(iR,1),N(1)));
iR_num=zeros(N(1),N(1),size(iR,3));
for ele=1:size(iR,3)
iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,ele),x_num,y_num);
end
clear x_ir y_ir x_num y_num
NF = [0*N; 0*N; 0*N];
nm=length(N);
numThetan=size(data_h,2);
angleRange=2;
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
Joint=-1; % 0: XRF; -1: XTM; 1: Joint inversion
ReconAttenu = 1; % 0: Recover W; 1: Recover miu
LogScale=1; %% determine if the XTM is solved taking log first or not
Beta=1e12; %% scaling parameter in joint inversion
%----------------------------------------------------------------------
% Compute the dependent-variable arrays
PlotObject=0;
plotSpec = 0; % Do you want to see the spectra? If so plotSpec = 1
plotTravel=0; % If plot the intersection of beam with object
plotUnit=0;
plotElement=0;
plotResult=0;
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
            XTM_Tensor;
        if(ReconAttenu)
            NumElement=1;
            W_level{level}=MU;
            if(level==1)
               WS= MU; 
            end
        else
            W_level{level}=W;
            if(level==1)
                WS= W; 
            end
        end
    else
        % XRF_XTM_Gaussian;
        XRF_XTM_Tensor;
        W_level{level}=W;
        if(level==1)
            WS=W;
        end
    end
    xtm_level{level}=DisR;
    L_level{level}=L;
    GI_level{level}=GlobalInd;
    nTau_level(level)= nTau;
    m_level(level,:)=m;
    SigmaT{level}=SigMa_XTM;
    %-----------------------------------------
    if(Joint~=-1)  
        xrf_level{level}=XRF;
        SI_level{level}=SelfInd;
        SigmaR{level}=SigMa_XRF;
    end
end
nTol=N(1)^2*NumElement;
%---------------------------------------------
% Specify initial guess for optimization.
rng('default');
W0=WS(:);
v0=rand(size(W0));%sum(reshape(x0,current_n,current_n,NumElement).*repmat(MUe,[current_n,current_n,1]),3);
xinitial=v0;
