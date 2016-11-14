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
synthetic=1;
if(synthetic)
    % sample='circle'; % one element mainly testing self-absorption 
    % sample = 'checkboard';
    % sample = 'Phantom';
    sample = 'fakeRod';
    NumElement=3;
else
    % sample='Seed';
    sample='Rod';
    NumElement=3;
end
N=100;% [33 17 9 5 3];% 17 9];%[129 65  9 5];%
angleScale=2; %1: half angle; 2: full angle
if(synthetic)
    numThetan=70; % number of scanning angles/projections
    DecomposedElement=0;
else
    Tol = 1e-2;
    omega=[-2 2 -2 2]*Tol; % units: cm
    if(strcmp(sample,'Seed'))
        load ./data/ApsDataExtract/DogaSeeds/DownSampleSeeds441_elements.mat
        load RawSpectra_seed.mat
        if(angleScale==1)
            cutInd=floor(size(data_H,3)/2);
            data_H=data_H(:,:,1:cutInd);
            spectra=spectra(:,:,1:cutInd);
        end
        data = permute(data_H,[2,3,1]);
        clear data_H
        slice =[9 13 14 17];%5:27;%      
        data_h=[];
        ang_rate=20;
        tau_rate=1;
        for ele=1:size(data,1)
            data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
        end
        if(ndims(data_h)==2)
            data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
        end
        truncChannel=0;% (DecomposedElement==0)*1;
        data_xrt=squeeze(data_h(3,:,:)); %transimission data
        s_a=0;
        I0=max(data_xrt(:));% reshape(data_h(2,:,:),size(data_h,2),size(data_h,3));% 
        % DisR=sparse(data_xrt');
        load DisR_removeStr
        data_xrf_decom=data_h(slice,:,:);
        numThetan=size(data_h,2);
        nTau=size(data_h,3)-1;
        N = size(data_h,3);
        data_xrf_raw=permute(spectra(1:tau_rate:end,:,1:ang_rate:end),[2 3 1]);
        data_xrf_raw=sparse(reshape(double(data_xrf_raw),[size(data_xrf_raw,1),size(data_xrf_raw,2)*size(data_xrf_raw,3)]));
        [x_ir,y_ir]=meshgrid(1:size(iR,1));
        [x_num,y_num]=meshgrid(linspace(1,size(iR,1),N(1)));
        iR_num=zeros(N(1),N(1),length(slice));
        for ele=1:size(iR_num,3)
            iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,slice(ele)),x_num,y_num);
        end
        save('tomopytest.mat','data_xrf_decom');
    elseif(strcmp(sample,'Rod'))
        load ~/Documents/Research/APSdata/GlassRod/2dSlice/Slice30
        load ~/Documents/Research/APSdata/GlassRod/2dSlice/tomoRec_3
        slice_tot = [3 4 25 30]; %GlassRod%
        data_h=[];
        ang_rate=1;
        tau_rate=9;
        for ele=1:size(data,1)
            data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
        end
        thetan = linspace(1,360,size(data,2));
        thetan_real = thetan(1:ang_rate:end); 
        if(ndims(data_h)==2)
            data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
        end
        numThetan=size(data_h,2);
        nTau=size(data_h,3)-1;
        N = size(data_h,3);
        truncChannel=0;%(DecomposedElement==0)*0;
        I0=reshape(data_h(43,:,:),size(data_h,2),size(data_h,3));
        data_sa=data_h(40,:,:);%% Scattering data: s_a;
        data_ds=data_h(38,:,:); %% Downstream Transmission 
        s_a=0;
        if(s_a)
            data_xrt=data_sa;
        else
            data_xrt=data_ds;
        end
        DisR=sparse(squeeze(sum(data_xrt(:,:,:),1))');
        data_xrf_decom=[];
        data_xrf_decom(1,:,:)=data_h(slice_tot(1),:,:)+data_h(slice_tot(2),:,:);
        data_xrf_decom(2,:,:)=data_h(slice_tot(3),:,:);
        data_xrf_decom(3,:,:)=data_h(slice_tot(4),:,:);
        % save('tomopytest.mat','data_xrf_decom','data_xrt');
        load spectra_30_aligned;
        data_xrf_raw=permute(spectra_30_aligned(1:tau_rate:end,1:ang_rate:end,:),[3 2 1]);
        data_xrf_raw=sparse(reshape(double(data_xrf_raw),[size(data_xrf_raw,1),size(data_xrf_raw,2)*size(data_xrf_raw,3)]));
        clear spectra_30_aligned data_sa data_ds data_h
        Si=0;
        W_element=0;
        if(W_element)
            data_xrf_decom=data_xrf_decom(3,:,:);
            iR=iR(:,:,3);
        end
        m_h=size(iR,1);
        [x_ir,y_ir]=meshgrid(1:m_h);
        [x_num,y_num]=meshgrid(linspace(1,m_h,N(1)));
        iR_num=zeros(N(1),N(1),size(iR,3));
        for ele=1:size(iR,3)
            iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,ele),x_num,y_num);
        end
    end
    clear data_xrt x_ir y_ir x_num y_num iR data spectra 
end
%%=============================
NoSelfAbsorption=0; % 0: include self-absorption in the XRF inversion
bounds = 1;  % no bound constraints
Joint=-1; % 0: XRF; -1: XTM; 1: Joint inversion
ReconAttenu = 1*(Joint==-1); % 0: Recover W; 1: Recover miu
Alternate=1*(Joint~=-1);
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
                W=MU;
            end
        else
            Forward_real;%XRF_XTM_Simplified;
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
                % [L,dd,~,~]=build_weight_matrix(MU,thetan,1,'area');
                XTM_Tensor;
                Lap{level}=delsq(numgrid('S',N(level)+2));
                L_level{level}=L;
            else
                IhH{level} = downdate_L(N, level);
                Lap{level}=IhH{level}*Lap{level-1}*IhH{level}';
                L_level{level}=L_level{level-1}*IhH{level}';
            end
            if(ReconAttenu)
                NumElement=1;
                if(level==1)
                    W_level{level}=MU;
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
        if(ReconAttenu)
            xtm_level_true{level}=DisR_true;
            xtm_level_pert{level}=DisR_pert;
        else
            xtm_level_true{level}=DisR;
        end
        GI_level{level}=GlobalInd;
        nTau_level(level)= nTau;
        m_level(level,:)=[N(level),N(level)];
        SigmaT{level}=SigMa_XTM;
        %-----------------------------------------
        if(Joint~=-1)
            if(synthetic)
                xrf_level{level}=XRF;
            else
                xrf_level_decom{level}=XRF_decom;
                xrf_level_raw{level}=XRF_raw;
            end
            SI_level{level}=SelfInd;
            SigmaR{level}=SigMa_XRF;
        end
    end
    clear L L_pert L_true W DisR_true DisR_pert GlobalInd SigMa_XTM SelfInd XRF XRF_decom XRF_raw SigMa_XRF data_h;
end;
nTol=N(1)^2*NumElement;
plotSim=0;
if(plotSim & Joint~=-1& n_level==1)
     a=squeeze(XRF_decom(:,:,2));b=squeeze(XRF_Simulated_decom(:,:,2));
     figure, for i=1:3; subplot(2,5,i);imagesc(XRF_Simulated_decom(:,:,i));colorbar; subplot(2,5,5+i);imagesc(XRF_decom(:,:,i));colorbar; end;
     a=squeeze(sum(sum(XRF_Simulated_raw,1),2));
     b=squeeze(sum(sum(XRF_decom{1},1),2)); 
     subplot(2,5,4),plot(DetChannel,squeeze(a),'r.-');
     ub=max(a(:));
     hold on;
     for ele=1:NumElement
         for el =1:4
         plot([BindingEnergy(ele,el) BindingEnergy(ele,el)],[0 ub],'g--');
         text(BindingEnergy(ele,el),ub,Element(Z(ele)));
         end
     end
     hold off;
     axis([0 12.1 0 ub+1e1]);
     subplot(2,5,9),plot(DetChannel,squeeze(b),'r.-');
     axis([0 12.1 0 max(b(:))+1e1]);
     subplot(2,5,5),semilogy(DetChannel(truncInd),a(truncInd),'r.-');
     axis([0 12.1 0 1e5]);
     subplot(2,5,10),
     semilogy(DetChannel(truncInd),b(truncInd),'r.-');
     axis([0 12.1 0 1e5]);
     err_xrf(scale+1)=norm(a(:)-b(:));
end

