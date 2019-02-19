global numThetan maxiter W0 sinoS
global totN xiter fiter ErrIter alpha

synthetic=0;
ReconAttenu = 0; % 0: Recover W; 1: Recover miu
initialize=1;
ind_scan=1;
slice=1;
do_setup;
rng('default');

ind_sub=4;
Z=Z(ind_sub);
NumElement=length(Z);
XRF_decom3D=XRF_decom3D(:,:,ind_sub,:);
maxiter=15;
totN=N^2*NumElement*nslice+numThetan*2;
Mt=permute(XRF_decom3D,[4 1 2 3]);
alpha=ones(nslice,NumElement);
Mt(isnan(Mt))=0;
for slice=1:nslice
    for ele=1:NumElement
        alpha(slice,ele)=std(reshape(Mt(:,:,slice,ele),(nTau+1)*numThetan,1));
    end
    alpha(slice,:)=alpha(slice,:)./sum(alpha(slice,:));
    % Mt(:,:,slice)=Mt(:,:,slice)./max(max(Mt(:,:,slice)));
end
for n=1:numThetan
    Mt(n,:,:)=Mt(n,:,:)./max(max(Mt(n,:,:)));
end
% Mt=Mt./max(Mt(:))+abs(min(Mt(:)));%%==normalize data;
% alpha=ones(size(alpha));

NF = [0*N; 0*N; 0*N];
Lmap=sparse(L);
MtPad=GaussianPadded(Mt);
    bounds=1;
    low=zeros(totN,1);
    up=inf*ones(totN,1);
    low(1:2:numThetan*2)=inf*floor(-nTau/2*ones(numThetan,1));
    up(1:2:numThetan*2)=inf*floor(nTau/2*ones(numThetan,1));
    low(2:2:numThetan*2)=inf*floor(-nslice/2*ones(numThetan,1));
    up(2:2:numThetan*2)=inf*floor(nslice/2*ones(numThetan,1));
    %%=================================================
    % fctn=@(x)sfun3D(x,Mt,Lmap);% on attenuation coefficients miu;
    % n_delta=2*numThetan;
    % W0=sparse(totN-n_delta,1);
    % err0=norm(W0-x0(n_delta+1:end));
    % [x,f,g]=tnbc(x0(n_delta+1:end),fctn,low(n_delta+1:end),up(n_delta+1:end));
    %%=================================================
    W0=sparse(totN,1);
    x0 = zeros(totN,1);
    % load manual_paunesku_shift_xy
    % x0(1:2:2*numThetan)=shift(:,1);
    % x0(2:2:2*numThetan)=shift(:,2);
    err0=norm(W0-x0);
    fctn_COR=@(x)sfun_shift_xy(x,MtPad,Lmap);% on attenuation coefficients miu;
    [f0,g0,aligned]=feval(fctn_COR,x0);
    if(bounds)
        [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
    else
        [xCOR,f,g,ierror] = tn (x0,fctn_COR);
    end
    [f,g,y]=feval(fctn_COR,xCOR);
    % save(['result/',sample,'/opt3D_xy_weighted',num2str(ang_rate),'_',num2str(tau_rate),char(Element(Z)),'noB.mat'],'xCOR','y');

