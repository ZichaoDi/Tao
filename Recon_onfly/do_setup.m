%---------------------------------------------------------------
% Set up a problem for optimization [multigrid or regular].
%---------------------------------------------------------------
global N numThetan NF
global bounds LogScale Joint
global grad_type err0 
global onlyXRF NoSelfAbsorption 
global W_level xrf_level_decom xrf_level_raw xtm_level L_level GI_level SI_level SigmaR SigmaT m_level nTau_level
global NumElement xinitial current_n frame I0 s_a
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
    % sample='circle'; % one element mainly testing self-absorption 
    % sample = 'checkboard';
    % sample = 'Phantom';
    sample = 'fakeRod';
    NumElement=3;
else
    % sample='Seed';
    sample='Xtal1';
    NumElement=3;
    % sample='Filter';
    % NumElement=6;
end
angleScale=1; %1: half angle; 2: full angle
if(synthetic)
    N=50;% [33 17 9 5 3];% 17 9];%[129 65  9 5];%
    numThetan=73; % number of scanning angles/projections
    DecomposedElement=0;
else
    Tol = 1e-2;
    omega=[-2 2 -2 2]*Tol; % units: cm
    if(strcmp(sample,'Seed'))
        setup_seed;
    elseif(strcmp(sample,'Filter'))
        setup_filter_sim;
    elseif(strcmp(sample,'Xtal1'))
        setup_Xtal1;
    elseif(strcmp(sample,'Rod'))
        setup_rod;
    end
    clear data_xrf_raw data_xrf_decom data_xrt x_ir y_ir x_num y_num iR data spectra 
end
%%=============================
NoSelfAbsorption=0; % 0: include self-absorption in the XRF inversion
bounds = 1;  % no bound constraints
Joint=1; % 0: XRF; -1: XTM; 1: Joint inversion
ReconAttenu = 1*(Joint==-1); % 0: Recover W; 1: Recover miu
Alternate=1*(Joint~=-1);
frame='LS';
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
onlyXRF=0;
%%------------------------------ Use same finest data for each level
%%-----------------------------------------------
current_n=N(1);
if(onlyXRF)
    SimulateXRF;
else
    if(Joint==-1)
        XTM_Tensor;
        if(ReconAttenu)
            NumElement=1;
            W=MU;
        end
    else
    end
end
nTol=N(1)^2*NumElement;
