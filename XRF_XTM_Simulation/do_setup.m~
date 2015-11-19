%---------------------------------------------------------------
% Set up a problem for optimization [multigrid or regular].
%---------------------------------------------------------------
global N numThetan NF
global bounds LogScale Joint
global grad_type err0 WS Beta W0
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
NoSelfAbsorption=0;
onlyXRF=0;
N=10;% [33 17 9 5 3];% 17 9];%[129 65  9 5];% 
NF = [0*N; 0*N; 0*N];
nm=length(N);
numThetan=5; % number of scanning angles/projections
numChannel=500;% number of energy channels on xrf detector
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
Joint=0; % 0: XRF; -1: XTM; 1: Joint inversion
Alternate=1;
LogScale=1; %% determine if the XTM is solved taking log first or not
% Beta=10^8;
%----------------------------------------------------------------------
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

if(level==1)
    % W=Phantom33;
else
W=downdate(W(:),1);
W=reshape(W,N(level),N(level),NumElement);
end
    current_n=N(level);
    XRF_XTM_Tensor;
%     XTM_Tensor;
    W_level{level}=W;
    if(level==1)
        WS=W;%MU;% 
    end
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
