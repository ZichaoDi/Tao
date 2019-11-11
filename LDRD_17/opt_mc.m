global N_delta maxiter W0 sinoS
global xiter fiter ErrIter
global alpha lambda
global reg_str

% ReconAttenu = 0; % 0: Recover W; 1: Recover miu
% slice=1;
% ind_cr=1;
% synthetic=0;
% sample='paunesku';
% do_setup;
% N_delta=numThetan;
% return;
rng('default');
res=3;
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,1);
initial_direction=size(d0,2);
deltaStar=zeros(N_delta,1);
W0=[deltaStar;W(:)];
maxiter=3;
aligned3D=zeros(numThetan*(nTau+1),NumElement,nslice);
recon=sparse(N^2*NumElement+N_delta,nslice);
for slice= 86;%1:nslice
    NF = [0*N; 0*N; 0*N];
    Mt=XRF_raw(:,:,:,slice);
    Mt=GaussianPadded(Mt,0,0);
    Mt(Mt<0)=0;
    prj=0*Mt;
    rng('default')
    snr=zeros(NumElement,1);
    for ele=1:NumElement
        poiss_ratio=2^alp;
        Mt(:,:,ele)=poissrnd(Mt(:,:,ele)*1/poiss_ratio);
        % gauss_ratio=7;
        % snr(ele)=sqrt(mean(mean(Mt(:,:,ele))));
        % Mt(:,:,ele)=Mt(:,:,ele)+randn(size(Mt(:,:,ele)))*gauss_ratio*1e-2*mean(mean(Mt(:,:,ele)));
        Mt(:,:,ele)=Mt(:,:,ele)+abs(min(min(Mt(:,:,ele))));
        prj(:,:,ele)=Mt(:,:,ele);
        a=Mt(:,:,ele);alpha(ele)=max(a(:));
        Mt(:,:,ele)=Mt(:,:,ele)./max(max(Mt(:,:,ele)));
    end
    % save(['pydataPoiss',num2str(poiss_ratio),'.mat'],'prj','thetan');
    alpha=alpha./sum(alpha);
    % alpha=ones(size(alpha));
    sinoS=reshape(Mt,numThetan*(nTau+1),NumElement);
    Lmap=sparse(L);
    for res_step=1:initial_direction
        delta=repmat(d0(:,res_step),n_delta/2,1);
        x0 = [(cos(theta)-1)*delta(1)+sin(theta)*delta(2);0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
        err0=norm(W0-x0);
        fctn_COR=@(x)sfun_shift_mc(x,Mt,Lmap);% on attenuation coefficients miu;
        % low=[-nTau/2*ones(N_delta,1);zeros(prod(m)*NumElement,1)];
        % up=[nTau/2*ones(N_delta,1);inf*ones(prod(m)*NumElement,1)];
        low=[-inf*ones(N_delta,1);zeros(prod(m)*NumElement,1)];
        up=[inf*ones(N_delta,1);inf*ones(prod(m)*NumElement,1)];

        bounds=1;
        % load aligned86_5_paunesku;
        % x0(1:numThetan)=shift;
        reg_str={'TV'};
        lam=-6;%[-15:2:-4];%, 0.2 0.4 0.8 1.6];
        for i=1:length(lam)
            lambda=0*10^lam(i);
            [~,~,alignedSignal]=feval(fctn_COR,x0);
            if(bounds)
                [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
                errW(res_step)=norm(W0(n_delta+1:end)-xCOR(n_delta+1:end));
                f_global(res_step)=f;
                [~,~,alignedSignal]=feval(fctn_COR,xCOR);
            else
                [xCOR,f,g,ierror] = tn (x0,fctn_COR);
            end
        end
    end
    recon(:,slice)=xCOR;
    aligned3D(:,:,slice)=alignedSignal;
end
prj=reshape(alignedSignal,numThetan,nTau+1,NumElement);
% save(['result/',sample,'/Recon_shift_mc_',num2str(numThetan),'_',num2str(NumElement),'_',num2str(nTau+1),'poiss',num2str(poiss_ratio),'.mat'],'xCOR','prj');

nrow=3;
figure,
xc=reshape(xCOR(N_delta+1:end,:),[N,N,NumElement,initial_direction]);
for ele=1:NumElement
    subplot(nrow,NumElement,ele);
    imagesc(Mt(:,:,ele));
    title(Element{Z(ele)});
set(gca,'xtick',[],'ytick',[]);
    if(ele==1);xlabel('\tau');ylabel('original: \theta');end
    subplot(nrow,NumElement,NumElement+ele);
    align=reshape(alignedSignal,numThetan*(nTau+1),NumElement);
    imagesc(reshape(align(:,ele),numThetan,nTau+1));
set(gca,'xtick',[],'ytick',[]);

    if(ele==1);xlabel('\tau');ylabel('aligned: \theta');end
    subplot(nrow,NumElement,NumElement*2+ele);
    imagesc(xc(:,:,ele));
    if(ele==1);ylabel('reconstruction');end
end

