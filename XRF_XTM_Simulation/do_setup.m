%---------------------------------------------------------------
% Set up a problem for optimization [multigrid or regular].
%---------------------------------------------------------------
% W,xrfData,xtmData,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0
%

global N numThetan NF
global bounds LogScale Joint
global grad_type err0 WS Beta W0
global NoSelfAbsorption
global W_level xrf_level xtm_level L_level GI_level SI_level SigmaR SigmaT m_level nTau_level
global DecomposedElement NumElement xinitial current_n
%--------------------------------------------------
% Select technique for gradient calculation.
%

grad_type = 'adj';  % 'adj' = adjoint/exact
% 'fdr' = finite-difference [real]
% 'fdi' = finite-difference [complex]

%--------------------------------------------
% Initialize arrays for discretizations
Tomo_startup;
%%-----------------------------------------------
load ./data/ApsDataExtract/DogaSeeds/DownSampleSeeds56_elements.mat
load tomoRecon56_half
iR=permute(out,[2 3 1]);
data=[];
slice=[8,9,14,17];
ang_rate=300;
for ele=1:length(slice)
data(:,ele,:)=sum(data_H(:,slice(ele),1:ang_rate:end),2);
end
data=data(:,:,1:floor(size(data,3)/2));
DecomposedElement=1;
% data=h5read('~/Documents/MATLAB/APSdata/xfm_Doga/xfm_data_elements.h5','/exchange/data');
% data=squeeze(sum(data,2));
% slice=1;
%%-----------------------------------------------
NoSelfAbsorption=0;
N=size(iR,1);%(floor(sqrt(size(data,1)^2/2));% 17 9];%[129 65  9 5];% 
NF = [0*N; 0*N; 0*N];
nm=length(N);
numThetan=size(data,3);
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
LogScale=1; %% determine if the XTM is solved taking log first or not
Beta=10^0;
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
else
    XRF_XTM_Tensor;
end
    W_level{level}=W;
    if(level==1)
        WS=  W;% MU;%
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
% NumElement=1;
nTol=N(1)^2*NumElement;
%---------------------------------------------
% Specify initial guess for optimization.
rng('default');
% load x_new;
W0=WS(:);
% x0=1e-1*rand(nTol/NumElement,1);%sum(reshape(x0,current_n,current_n,NumElement).*repmat(MUe,[current_n,current_n,1]),3);
%{
xinitial=x0;
x0 = 10^(-1)*rand(nTol,1);%reshape(iR(:,:,slice),[size(iR,1)*size(iR,2),1]);%10^(-1)*rand(nTol,1)+WS(:);%ww(:);%;%zeros(size(WS(:)));%
err0=norm(x0-WS(:));
v0=x0;
%}
