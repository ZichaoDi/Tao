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
% n_delta=2*1;%numThetan;
%%===============Load Sample=====================
% synthetic=0;
if(synthetic)
    % sample='circle'; % one element mainly testing self-absorption 
    % sample='Golosio';
    % sample = 'checkboard';
    % sample = 'MRI';
    % sample = 'Phantom';
    % sample = 'fakeRod';
    NumElement=1;
else
    % sample='Run02';%'olga';%'miller';%'Rod';%'Filter';%
    % sample='Paunesku';%'miller';%'Zn_modified_6';
    sample='chip';
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
    % numThetan=30; % number of scanning angles/projections
    DecomposedElement=0;
else
    Tol = 1e-2;
    omega=[-2 2 -2 2]*Tol; % units: cm
    if(strcmp(sample,'Seed'))
        setup_seed;
    elseif(strcmp(sample,'synthetic'))
        load syntheticMultiSino.mat; 
        thetan_real=thetan;
        XRF_raw_tot=prj;
        XRF_raw=prj;
        N=size(wtrue,1);
        numThetan=size(prj,1);
        nTau=size(prj,2)-1;
        nslice=1;
        Ztot=[20 25 36];
        % ele_ind=1;
        Z=Ztot(ele_ind);
    elseif(strcmp(sample,'Filter'))
        setup_filter;
    elseif(strcmp(sample,'Rod'))
        setup_rod;
    elseif(strcmp(sample,'Paunesku'))
        setup_paunesku;
    elseif(strcmp(sample,'miller'))
        setup_miller;
    elseif(strcmp(sample,'20173'))
        setup_20173;
    elseif(strcmp(sample,'Run02'))
        nslice=1200;
        load(['Run02_',num2str(slice)]);
        XRF_decom=-log(data(1:tau_rate:end,1:ang_rate:end)'./max(max(data)));
        nTau=size(XRF_decom,2)-1;
        N=nTau+1;
        numThetan=size(XRF_decom,1);
        thetan_real=linspace(0,180,numThetan);
    elseif(strcmp(sample,'olga'))
        data=readNPY('XRF data/olga-dat.npy');
        ang=readNPY('XRF data/olga-ang.npy');
        ang_range=[1:ang_rate:length(ang)];
        [thetan_real,ind_ang]=sort(double(ang(1:ang_rate:end)'));
        numThetan=length(thetan_real);
        nslice=size(data,3);
        XRF_decom3D=double(permute(squeeze(data([30 37],ang_range(ind_ang),:,1:tau_rate:end)),[2 3 4 1]));
        nTau=size(XRF_decom3D,3)-1;
        N=nTau+1;
        Z=[26 34];
    elseif(strcmp(sample,'Zn_modified_6'))
        data=[];
        numThetan=22;
        files=dir('~/Documents/Research/APS/Zn_modified_6/*.tif');
        for n=1:numThetan
            [pathstr,name,ext]=fileparts(['~/Documents/Research/APS/Zn_modified_6/',files(n).name]);
            data(n,:,:)=imread([pathstr,'/', name, ext]);
        end
        XRF_decom3D=data;%permute(data,[1 3 2]);
        nslice=size(XRF_decom3D,2);
        nTau=size(XRF_decom3D,3)-1;
        N=nTau+1;
        thetan_real=importdata('angle.txt');
        thetan_real=thetan_real';
    else
        run(['setup_',sample,'.m']);
    end
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
frame='EM';
%%--------------------------------------------------------------------

XTM_Tensor;
if(ReconAttenu)
    NumElement=1;
    W=MU_XTM;
end
%-----------------------------------------
nTol=N(1)^2*NumElement;
