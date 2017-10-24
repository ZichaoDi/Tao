%---------------------------------------------------------------
% Set up a problem for optimization [multigrid or regular].
%---------------------------------------------------------------
global N numThetan NF nslice
global bounds LogScale Joint
global grad_type err0 
global onlyXRF NoSelfAbsorption coarsen_type 
global level W_level xrf_level_decom xrf_level_raw xtm_level L_level GI_level SI_level SigmaR SigmaT m_level nTau_level
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
    % sample='Golosio';
    % sample = 'checkboard';
    sample = 'Phantom';
    % sample = 'fakeRod';
    NumElement=6;
else
    sample='Paunesku';%'miller';%'Rod';%'Filter';%
    NumElement=6;
end
coarsen_type='smooth';
N=[65];%[33 17 9 5 3];% 17 9];%[129 65  9 5];%
angleScale=2; %1: half angle; 2: full angle
if(synthetic)
    numThetan=10; % number of scanning angles/projections
    DecomposedElement=0;
else
    Tol = 1e-2;
    omega=[-2 2 -2 2]*Tol; % units: cm
    if(strcmp(sample,'Seed'))
        setup_seed;
    elseif(strcmp(sample,'Filter'))
        setup_filter;
    elseif(strcmp(sample,'Rod'))
        setup_rod;
    elseif(strcmp(sample,'Paunesku'))
        setup_paunesku;
    elseif(strcmp(sample,'miller'))
        setup_miller;
    end
    clear data_xrf_raw data_xrf_decom data_xrt x_ir y_ir x_num y_num iR data spectra 
end
%%=============================
NoSelfAbsorption=0; % 0: include self-absorption in the XRF inversion
bounds = 1;  % no bound constraints
Joint=-1; % 0: XRF; -1: XTM; 1: Joint inversion
Alternate=1*(Joint~=-1);
% ReconAttenu = 0*(Joint==-1); % 0: Recover W; 1: Recover miu
frame='EM';
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
n_level=length(N);
if(n_level==1)
    current_n=N(1);
    if(onlyXRF)
        SimulateXRF;
    else
        if(Joint==-1)
            XTM_Tensor;
            if(ReconAttenu)
                NumElement=1;
                W=MU_XTM;
            end
            % L=build_weight_matrix(W,thetan,0);
            % L=reshape(permute(reshape(L,nTau+1,numThetan,N^2),[2 1 3]),(nTau+1)*numThetan,N^2);
        else
            Forward_real;% XRF_XTM_Simplified; %
        end
    end
    %-----------------------------------------
    clear L_pert L_true DisR_true DisR_pert;
else
    nh=N(1);
    level=[1:n_level];
    N=nh./(2.^(level-1))+(2.^(level-1)-1)./2.^(level-1);
    NF = [0*N; 0*N; 0*N];
    nm=length(N);
    W_level=cell(nm,1);
    if(synthetic)
        xrf_level=cell(nm,1);
    else
        xrf_level_decom=cell(nm,1);
        xrf_level_raw=cell(nm,1);
    end
    xtm_level=cell(nm,1);
    L_level=cell(nm,1);
    GI_level=cell(nm,1);
    SI_level=cell(nm,1);
    SigmaR=cell(nm,1);
    SigmaT=cell(nm,1);
    m_level=zeros(nm,2);
    nTau_level= zeros(nm,1);
    for level=1:nm
        current_n=N(level);
        if(Joint==-1)
            if(level==1)
                XTM_Tensor;
                L_level{level}=L;%build_weight_matrix(W,thetan,0);
            else
                L_level{level}=downdate_radon(L_level{level-1},numThetan,nTau);
            end
            if(ReconAttenu)
                NumElement=1;
                if(level==1)
                    W_level{level}=MU_XTM;
                else
                    W_level{level}=N(level)^2;
                end
            else
                W_level{level}=W;
            end
        else
            XRF_XTM_Simplified;
            W_level{level}=W;
            L_level{level}=L;
        end
        xtm_level_true{level}=DisR_Simulated;
        nTau_level(level)= nTau;
        m_level(level,:)=[N(level),N(level)];
        %-----------------------------------------
        if(Joint~=-1)
            if(synthetic)
                xrf_level{level}=XRF;
            else
                xrf_level_decom{level}=XRF_decom;
                xrf_level_raw{level}=XRF_raw;
            end
        end
    end
    clear L L_pert L_true W DisR_true DisR_pert GlobalInd SigMa_XTM SelfInd XRF XRF_decom XRF_raw SigMa_XRF data_h;
end;
nTol=N(1)^2*NumElement;
