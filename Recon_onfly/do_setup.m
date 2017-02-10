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
    sample='Rod';
    NumElement=3;
end
angleScale=2; %1: half angle; 2: full angle
if(synthetic)
    N=50;% [33 17 9 5 3];% 17 9];%[129 65  9 5];%
    numThetan=73; % number of scanning angles/projections
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
        ang_rate=2;
        tau_rate=20;
        for ele=1:size(data,1)
            data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
        end
        if(ndims(data_h)==2)
            data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
        end
        data_xrt=squeeze(data_h(3,:,:)); %transimission data
        % I0=max(data_xrt(:));
        I0=reshape(data_h(2,:,:),size(data_h,2),size(data_h,3));% 
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
    elseif(strcmp(sample,'Filter'))
        Filter_setup;
    elseif(strcmp(sample,'Rod'))
        % load ~/Documents/Research/APS/GlassRod/2dSlice/Slice30
        % ind_i0=43;
        % ind_xrt=38;
        % slice_tot = [3 4 25 30]; %GlassRod%
        load hong30;
        data=MAPS;
        ind_i0=33;
        ind_xrt=32;
        slice_tot = [2 3 16 20];
        %%%================================
        c_sh=floor(size(data,3)/2-842);
        data_temp=data;
        center_shift=c_sh;
        if(c_sh>=0)
            for sub_ch=1:size(data,1)
                data_temp(sub_ch,:,center_shift:end)=data(sub_ch,:,1:end-center_shift+1);
                data_temp(sub_ch,:,1:center_shift-1)=data(sub_ch,:,end-center_shift+2:end);
            end
        else
            for sub_ch=1:size(data,1)
                data_temp(sub_ch,:,1:end+center_shift)=data(sub_ch,:,abs(center_shift)+1:end);
                data_temp(sub_ch,:,end+center_shift+1:end)=data(sub_ch,:,1:abs(center_shift));
            end
        end
        data=data_temp;
        clear data_temp;
        %%===========================filter original silicon detector XRT
        windowSize=10;
        xrt1=zeros(size(data(ind_xrt,:,:)));
        for i=1:73, 
            % y=filter(1/windowSize*ones(1,windowSize),1,[squeeze(data_xrt(1,i,1))*ones(windowSize,1);squeeze(data_xrt(1,i,:))]); 
            % xrt1(1,i,:)=y(11:end);end
            if(i==72 | i==73)
                xrt1(1,i,:)=medfilt1(data(ind_xrt,i,:),windowSize,'truncate')-30;
            else
                xrt1(1,i,:)=medfilt1(data(ind_xrt,i,:),windowSize,'truncate');
            end
        end
        data(ind_xrt,:,:)=xrt1;
        clear y xrt1;

        load tomoRod
        data_h=[];
        ang_rate=2;
        tau_rate=20;
        for ele=1:size(data,1)
            data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
            % data_h(ele,:,:)=sum_interval(squeeze(data(ele,:,:)),ang_rate,tau_rate);
        end
        thetan = linspace(1,360,size(data,2));
        thetan_real = thetan(1:ang_rate:end); 
        if(ndims(data_h)==2)
            data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
        end
        numThetan=size(data_h,2);
        thetan=linspace(0,360,numThetan);
        nTau=size(data_h,3)-1;
        N = size(data_h,3);
        data_xrt=data_h(ind_xrt,:,:); %% Downstream Transmission 
        %%==============================================
        I0=max(data_xrt(:));
        % I0=reshape(data_h(ind_i0,:,:),size(data_h,2),size(data_h,3));
        DisR=sparse(squeeze(sum(data_xrt(:,:,:),1))');
        data_xrf_decom=[];
        data_xrf_decom(1,:,:)=data_h(slice_tot(1),:,:)+data_h(slice_tot(2),:,:);
        data_xrf_decom(2,:,:)=data_h(slice_tot(3),:,:);
        data_xrf_decom(3,:,:)=data_h(slice_tot(4),:,:);
        % save('tomopy_coarse.mat','data_xrf_decom','data_xrt');
        XRF_decom=permute(data_xrf_decom(:,:,:),[2 3 1]);
        %%================================================
        % load spectra_30_aligned;
        load hong_spectra30;
        spectra_30_aligned=spectra;
        %%================================================
        spectra=double(0.*spectra_30_aligned);
        center_shift=c_sh;
        if(center_shift>=0)
            spectra(center_shift:end,:,:)=spectra_30_aligned(1:end-center_shift+1,:,:);
            spectra(1:center_shift-1,:,:)=spectra_30_aligned(end-center_shift+2:end,:,:);
        else
            spectra(1:end+center_shift,:,:)=spectra_30_aligned(abs(center_shift)+1:end,:,:);
            spectra(end+center_shift+1:end,:,:)=spectra_30_aligned(1:abs(center_shift),:,:);
        end

        % data_xrf_raw=permute(sum_interval(spectra,tau_rate,ang_rate),[3 2 1]);
        data_xrf_raw=permute(spectra(1:tau_rate:end,1:ang_rate:end,:),[3 2 1]);
        data_xrf_raw=sparse(reshape(double(data_xrf_raw),[size(data_xrf_raw,1),size(data_xrf_raw,2)*size(data_xrf_raw,3)]));
        XRF_raw=data_xrf_raw';
        %%====== without shifting
        % data_xrf_raw=permute(spectra_30_aligned(1:tau_rate:end,1:ang_rate:end,:),[3 2 1]);
        % data_xrf_raw=sparse(reshape(double(data_xrf_raw),[size(data_xrf_raw,1),size(data_xrf_raw,2)*size(data_xrf_raw,3)]));
        % XRF_raw=data_xrf_raw';
        %%=================================
        clear spectra spectra_30_aligned data_sa data_ds data_h
        m_h=size(iR,1);
        [x_ir,y_ir]=meshgrid(1:m_h);
        [x_num,y_num]=meshgrid(linspace(1,m_h,N(1)));
        iR_num=zeros(N(1),N(1),size(iR,3));
        for ele=1:size(iR,3)
            iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,ele),x_num,y_num);
        end
    end
    clear data_xrf_raw data_xrf_decom data_xrt x_ir y_ir x_num y_num iR data spectra 
end
%%=============================
NoSelfAbsorption=0; % 0: include self-absorption in the XRF inversion
bounds = 1;  % no bound constraints
Joint=1; % 0: XRF; -1: XTM; 1: Joint inversion
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
