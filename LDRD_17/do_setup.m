%---------------------------------------------------------------
% Set up a problem for optimization [multigrid or regular].
%---------------------------------------------------------------
global N numThetan NF nslice
global bounds LogScale Joint
global nchannel grad_type err0 
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
% n_delta=2*1;%numThetan;
%%===============Load Sample=====================
% synthetic=1;
if(synthetic)
    % sample='circle'; % one element mainly testing self-absorption 
    % sample='Golosio';
    % sample = 'checkboard';
    % sample = 'MRI';
    % sample = 'Phantom';
    % sample = 'fakeRod';
    % N=50;
    % numThetan=30; % number of scanning angles/projections
    NumElement=1;
else
    % sample='run02';%'olga';%'miller';%'rod';%'filter';%
    % sample='paunesku';%'miller';%'Zn_modified_6';
    % sample='chip';
    % sample = 'synthetic';
end
if(~exist('slice','var'))
    slice=1;
end
if(~exist('element_ind','var'))
    element_ind=1;
end
initialize=1;

coarsen_type='smooth';
angleScale=2; %1: half angle; 2: full angle
if(synthetic)
    DecomposedElement=0;
    if(strcmp(sample,'Golosio'))
            CreateCircle; 
            Z = 19;%[6 8 14 20 26];% Golosio's Sample
        elseif(strcmp(sample,'checkboard'))
            W = kron(invhilb(round(N(1)/10))<0, ones(10,10));
            Z = 19;
            W=repmat(W,[1 1 NumElement]);
    elseif(strcmp(sample,'Phantom')|strcmp(sample,'MRI'))
            Z = [19];
            CreateElement; 
    elseif(strcmp(sample,'circle'))
            [X,Y]=meshgrid(1:N,1:N);
            center=[N/3, N/3];
            r=5;
            pix = (X-center(1)).^2+(Y-center(2)).^2 <= r^2;%& (X-center(1)).^2+(Y-center(2)).^2>=(r-5).^2; %% circle
            W(pix)=1;
            Z = [14];
    elseif(strcmp(sample, 'fakeRod'))
            CreateRod;
            Z=[79 5 8 11 13 14 19 74];% Glass Rod
    end
else
    Tol = 1e-2;
    omega=[-2 2 -2 2]*Tol; % units: cm
    run(['setup_',sample,'.m']);
    W=zeros(N(1),N(1),NumElement);
    NumElement=length(Z);
    clear data_xrf_raw data_xrf_decom data_xrt x_ir y_ir x_num y_num iR data spectra 
    n_delta=2*numThetan;
    ind_scan=0;
end
%%=============================
bounds = 1;  % no bound constraints
Joint=-1; % 0: XRF; -1: XTM; 1: Joint inversion
if(~exist('ReconAttenu','var'))
    ReconAttenu = 1;% 0: Recover W; 1: Recover miu
end
frame='LS';
%%--------------------------------------------------------------------

XTM_Tensor;
if(ReconAttenu)
    NumElement=1;
    W=MU_XTM;
end
%-----------------------------------------
nTol=N(1)^2*NumElement;
