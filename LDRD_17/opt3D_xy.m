global numThetan maxiter W0 sinoS
global totN xiter fiter ErrIter alpha

ReconAttenu = 0; % 0: Recover W; 1: Recover miu
initialize=1;
ind_scan=1;
slice=1;
tau_rate=10;
ang_rate=10;
do_setup;
rng('default');

maxiter=100;
totN=N^2*NumElement*nslice+numThetan*2;
Mt=permute(XRF_decom3D,[3 2 4 1]);
alpha=ones(nslice,NumElement);
Mt(isnan(Mt))=0;
for slice=1:nslice
    for ele=1:NumElement
        alpha(slice,ele)=std(reshape(Mt(:,slice,ele,:),(nTau+1)*numThetan,1));
    end
    alpha(slice,:)=alpha(slice,:)./sum(alpha(slice,:));
end
Mt=Mt./max(Mt(:))+abs(min(Mt(:)));%%==normalize data;
alpha=ones(size(alpha));

NF = [0*N; 0*N; 0*N];
Lmap=sparse(L);
MtPad=GaussianPadded(Mt);
    x0 = zeros(totN,1);
    bounds=1;
    low=zeros(totN,1);
    up=inf*ones(totN,1);
    low(1:2:numThetan*2)=floor(-nTau/2*ones(numThetan,1));
    up(1:2:numThetan*2)=floor(nTau/2*ones(numThetan,1));
    low(2:2:numThetan*2)=floor(-nslice/2*ones(numThetan,1));
    up(2:2:numThetan*2)=floor(nslice/2*ones(numThetan,1));
    %%=================================================
    % load shift_manual
    % x0(1:2:n_delta)=-align(ang_ind,1);
    % x0(2:2:n_delta)=-align(ang_ind,2);
    % [~,~,aligned]=feval(fctn_COR,x0);
    fctn=@(x)sfun3D(x,permute(Mt,[4 1 2 3]),Lmap);% on attenuation coefficients miu;
    n_delta=2*numThetan;
    W0=sparse(totN-n_delta,1);
    err0=norm(W0-x0(n_delta+1:end));
    [x,f,g]=tnbc(x0(n_delta+1:end),fctn,low(n_delta+1:end),up(n_delta+1:end));
    %%=================================================
    W0=sparse(totN,1);
    x0=zeros(totN,1);
    err0=norm(W0-x0);
    fctn_COR=@(x)sfun_shift_xy(x,MtPad,Lmap);% on attenuation coefficients miu;
    if(bounds)
        [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
    else
        [xCOR,f,g,ierror] = tn (x0,fctn_COR);
    end
    % save(['result/',sample,'/opt3D_xy',num2str(ang_rate),'_',num2str(tau_rate),'.mat'],'xCOR','x');
