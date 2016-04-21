%---------------------------------------------------------------
% Set up a problem for optimization [multigrid or regular].
%---------------------------------------------------------------
global N numThetan NF
global bounds LogScale Joint
global grad_type err0 
global onlyXRF NoSelfAbsorption pix_inner
global W_level xrf_level_decom xrf_level_raw xtm_level L_level GI_level SI_level SigmaR SigmaT m_level nTau_level
global NumElement xinitial current_n EM
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
    % sample='circle'; % one element mainly testing self-absorption; 
    sample = 'checkboard';
    % sample = 'Phantom';
else
    % sample='Seed';
    sample='Rod';
end
N=100;% [33 17 9 5 3];% 17 9];%[129 65  9 5];%
angleScale=2; %1: half angle; 2: full angle
if(synthetic)
    numThetan=70; % number of scanning angles/projections
    numChannel=5;% number of energy channels on xrf detector
    DecomposedElement=0;
else
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
        slice =[13 14 17];%5:27;%      
        data_h=[];
        ang_rate=100;
        tau_rate=40;
        for ele=1:size(data,1)
            data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
        end
        if(ndims(data_h)==2)
            data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
        end
        DecomposedElement=1;
        truncChannel=1;% (DecomposedElement==0)*1;
        data_xrt=data_h(3,:,:); %transimission data
        data_xrf_decom=data_h(slice,:,:);
        data_xrf_raw=permute(spectra(1:tau_rate:end,:,1:ang_rate:end),[2 3 1]);
        [x_ir,y_ir]=meshgrid(1:size(iR,1));
        [x_num,y_num]=meshgrid(linspace(1,size(iR,1),N(1)));
        iR_num=zeros(N(1),N(1),size(iR,3));
        for ele=1:size(iR,3)
            iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,ele),x_num,y_num);
        end
        save('tomopytest.mat','data_xrf_decom');
    elseif(strcmp(sample,'Rod'))
        load ~/Documents/Research/APSdata/GlassRod/2dSlice/Slice30
        load ~/Documents/Research/APSdata/GlassRod/2dSlice/tomoRec_3
        slice_tot = [3 4 25 30]; %GlassRod%
        slice = [4 25 30];
        data_h=[];
        ang_rate=1;
        tau_rate=15;
        for ele=1:size(data,1)
            data_h(ele,:,:)=sum(data(ele,1:ang_rate:end,1:tau_rate:end),1);
        end
        thetan = linspace(1,360,size(data,2));
        thetan_real = thetan(1:ang_rate:end); 
        if(ndims(data_h)==2)
            data_h=reshape(data_h,size(data_h,1),1,size(data_h,2));
        end
        truncChannel=0;%(DecomposedElement==0)*0;
        I0=squeeze(data_h(43,:,:));
        data_xrt=data_h(38,:,:); 
        data_xrf_decom=[];
        data_xrf_decom(1,:,:)=data_h(slice_tot(1),:,:)+data_h(slice_tot(2),:,:);
        data_xrf_decom(2,:,:)=data_h(slice_tot(3),:,:);
        data_xrf_decom(3,:,:)=data_h(slice_tot(4),:,:);
        save('tomopytest.mat','data_xrf_decom','data_xrt');
        load spectra_30;
        data_xrf_raw=permute(spectra_30(1:tau_rate:end,1:ang_rate:end,:),[3 2 1]);
        clear spectra_30
        Tol = 1e-2;
        omega=[-2 2 -2 2]*Tol*2; % units: cm
        m_h=size(iR,1);
        xc=getCellCenteredGrid(omega,[m_h,m_h]);
        xc=reshape(xc,m_h^2,2);
        center = [1.4 -1.2 ; 0 0; 1.1 1.1].*Tol;
        center = [2.9 -2.7 ; 0 0; 2.6 2.6].*Tol;
        r = [ 1.7e-1 3.4 1.7e-1].*Tol;
        w_ele=[19.3 2.33 19.3]; % units: g/cm^3
        w_ele=w_ele.*[1 1 1]*2;
        iR_fill = zeros(m_h^2,3);
        for ele=1:3
            pix = (xc(:,1)-center(ele,1)).^2+(xc(:,2)-center(ele,2)).^2 <= r(ele)^2;
            if(ele==1)
                % iR_fill(find(pix==1),3)=w_ele(3);%max(reshape(iR(:,:,3),m_h^2,1));
                iR_fill(find(pix==1),1)=w_ele(1);%max(reshape(iR(:,:,1),m_h^2,1));
                % iR_fill(find(pix==1),2)=w_ele(2)+1;%6e-3;
            elseif(ele==2)
                % rec=[[1.1e-5 1.1e-5+r(3)];[1.2e-5 -1.2e-5+r(3)];[1.1e-5 1.1e-5-r(3)];[1.2e-5 -1.2e-5-r(3)]];
                % rec=[center(1,:);center(1,:)-[r(3)/2,0];center(3,:)-[r(3)/2,0];center(3,:)];
                % pix_extra=isPointInPolygon(xc,rec); 
                % pix = (xc(:,1)-center(ele,1)).^2+(xc(:,2)-center(ele,2)).^2 <= r(ele)^2 & (xc(:,1)-center(ele,1)).^2+(xc(:,2)-center(ele,2)).^2 >= (2*r(ele)/3)^2;
                iR_fill(find(pix==1),2)=w_ele(2);%6e-3;
                % iR_fill(find(pix_extra==1),2)=w_ele(2);%6e-3;
            elseif(ele==3)
            iR_fill(find(pix==1),ele)=w_ele(ele);%max(reshape(iR(:,:,ele),m_h^2,1));
            % iR_fill(find(pix==1),2)=w_ele(2)+1;%6e-3;%max(reshape(iR(:,:,2),m_h^2,1));
            end
        end
        % iR=reshape(iR_fill*1e0,[m_h,m_h,3]);
        Si=0;
        W_element=0;
        if(W_element)
            data_xrf_decom=data_xrf_decom(3,:,:);
            iR=iR(:,:,3);
        end
        iR= iR*1;%*10^(scale);
        clear iR_2;
        [x_ir,y_ir]=meshgrid(1:m_h);
        [x_num,y_num]=meshgrid(linspace(1,m_h,N(1)));
        iR_num=zeros(N(1),N(1),size(iR,3));
        for ele=1:size(iR,3)
            iR_num(:,:,ele)=interp2(x_ir,y_ir,iR(:,:,ele),x_num,y_num);
        end
    end
    clear x_ir y_ir x_num y_num iR data spectra 
    numThetan=size(data_h,2);
end
%%-----------------------------------------------
nh=N(1);
n_level=1;
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
%%=============================
NoSelfAbsorption=0; % 0: include self-absorption in the XRF inversion
bounds = 1;  % no bound constraints
Joint=1; % 0: XRF; -1: XTM; 1: Joint inversion
ReconAttenu = 0*(Joint==-1); % 0: Recover W; 1: Recover miu
Alternate=1*(Joint~=-1);
EM=1;
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
        XRF_XTM_Simplified;%Tensor;
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

nTol=N(1)^2*NumElement;
if(Joint~=-1)
     xrfs=XRF_Simulated_decom;xrfr=xrf_level_decom{1}; 
     a=squeeze(xrfr(:,:,2));b=squeeze(xrfs(:,:,2));
     figure, for i=1:3; subplot(2,5,i);imagesc(xrfs(:,:,i));colorbar; subplot(2,5,5+i);imagesc(xrfr(:,:,i));colorbar; end;
     a=squeeze(sum(sum(XRF_Simulated_raw,1),2));
     b=squeeze(sum(sum(xrf_level_raw{1},1),2)); 
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

% save([sample,num2str(N(1)),'_',num2str(numThetan),'_',num2str(nTau+1),'_',num2str(numChannel),'.mat'],'-v7.3')
