%---------------------------------------------------------------
% Set up a problem for optimization [multigrid or regular].
%---------------------------------------------------------------
% W,xrfData,xtmData,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0
%

global N
global bounds
global grad_type

%--------------------------------------------------
% Select technique for gradient calculation.
%

grad_type = 'adj';  % 'adj' = adjoint/exact
% 'fdr' = finite-difference [real]
% 'fdi' = finite-difference [complex]

%--------------------------------------------
% Initialize arrays for discretizations

N=[17];% 9 5 3];
nm=length(N);
W_level=cell(nm,1);
xrf_level=cell(nm,1);
xtm_level=cell(nm,1);
L_level=cell(nm,1);
GI_level=cell(nm,1);
SI_level=cell(nm,1);
m_level=zeros(nm,2);
nTau_level= zeros(nm,1);
bounds   = 1;  % no bound constraints
%----------------------------------------------------------------------
% Compute the dependent-variable arrays
PlotObject=0;
plotSpec = 0; % Do you want to see the spectra? If so plotSpec = 1
plotTravel=0; % If plot the intersection of beam with object
plotUnit=0;
plotElement=0;
% plotResult=0;
LogScale=1;
for level=1:nm
    current_n=N(level);
    XRF_XTM_Tensor;
    W_level{level}=W;
    xrf_level{level}=XRF;
    xtm_level{level}=DisR;
    L_level{level}=L;
    GI_level{level}=GlobalInd;
    SI_level{level}=SelfInd;
    m_level(level,:)=m;
    nTau_level(level)= nTau;
end
nTol=N(1)*NumElement;
%---------------------------------------------
% Specify initial guess for optimization.
rng('default');
x0=1*10^(-1)*rand(nTol,1);