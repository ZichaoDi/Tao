%---------------------------------------------------------------
% Set up a problem for optimization [multigrid or regular].
%---------------------------------------------------------------
% W,xrfData,xtmData,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0
%

global N numThetan
global bounds LogScale Joint
global grad_type err0 WS
global onlyXRF

%--------------------------------------------------
% Select technique for gradient calculation.
%

grad_type = 'adj';  % 'adj' = adjoint/exact
% 'fdr' = finite-difference [real]
% 'fdi' = finite-difference [complex]

%--------------------------------------------
% Initialize arrays for discretizations
Tomo_startup;
onlyXRF=0;
N=50;%[17 9 5 3];%[3];%
nm=length(N);
numThetan=10;%[2 2 1 1];
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
bounds   = 1;  % no bound constraints
Joint=1; % 0: XRF; -1: XTM; 1: Joint inversion
%----------------------------------------------------------------------
% Compute the dependent-variable arrays
PlotObject=0;
plotSpec = 0; % Do you want to see the spectra? If so plotSpec = 1
plotTravel=0; % If plot the intersection of beam with object
plotUnit=0;
plotElement=0;
% plotResult=1;
LogScale=1;
for level=1:nm
    current_n=N(level);
    XRF_XTM_Tensor;
    W_level{level}=W;
    if(level==1)
        WS=W;
    end
    xrf_level{level}=XRF;
    xtm_level{level}=DisR;
    L_level{level}=L;
    GI_level{level}=GlobalInd;
    SI_level{level}=SelfInd;
    m_level(level,:)=m;
    nTau_level(level)= nTau;
    SigmaT{level}=SigMa_XTM;
    SigmaR{level}=SigMa_XRF;
end
nTol=N(1)^2*NumElement;
%---------------------------------------------
% Specify initial guess for optimization.
rng('default');
x0=WS(:)+1*10^(-1)*rand(nTol,1);
xinitial=x0;
err0=norm(x0-WS(:));