global nchannel N_delta maxiter W0 sinoS
global xiter fiter ErrIter
global alpha lambda
global reg_str

% ReconAttenu = 0; % 0: Recover W; 1: Recover miu
% slice=1;
% synthetic=0;
% ele_ind=1;
% do_setup;
% N_delta=numThetan;
rng('default');
res=3;
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,1);
initial_direction=size(d0,2);
deltaStar=zeros(N_delta,1);
maxiter=300;
xtot=zeros(N^2+N_delta,nchannel);
for ele=1:nchannel;
    W0=[deltaStar;reshape(W(:,:,ele),N^2,1)];
    for slice= 86;%1:nslice
        NF = [0*N; 0*N; 0*N];
        MtS=Mt(:,:,ele);
        for res_step=1:initial_direction
            delta=repmat(d0(:,res_step),n_delta/2,1);
            x0 = [(cos(theta)-1)*delta(1)+sin(theta)*delta(2);0*10^(0)*rand(m(1)*m(2),1)];
            err0=norm(W0-x0);
            alpha=1e0;
            fctn_COR=@(x)sfun_shift_mc(x,MtS,Lmap);% on attenuation coefficients miu;
            low=[-nTau/2*ones(N_delta,1);zeros(prod(m),1)];
            up=[nTau/2*ones(N_delta,1);inf*ones(prod(m),1)];
            % low=[-inf*ones(N_delta,1);zeros(prod(m),1)];
            % up=[inf*ones(N_delta,1);inf*ones(prod(m),1)];
            bounds=1;
            % load aligned86_5_paunesku;
            % x0(1:numThetan)=shift;
            reg_str={'TV'};
            lam=-6;%[-15:2:-4];%, 0.2 0.4 0.8 1.6];
            for i=1:length(lam)
                lambda=0*10^lam(i);
                if(bounds)
                    [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
                    errW(res_step)=norm(W0(n_delta+1:end)-xCOR(n_delta+1:end));
                    f_global(res_step)=f;
                    [~,~,alignedSignal]=feval(fctn_COR,xCOR);
                    %% ============================================
                else
                    [xCOR,f,g,ierror] = tn (x0,fctn_COR);
                end
            end
        end
    end
    prj(:,:,ele)=reshape(alignedSignal,numThetan,nTau+1);
    xtot(:,ele)=xCOR;
end
nrow=3;
figure,
xc=reshape(xtot(N_delta+1:end,:),[N,N,ele]);
for i=1:ele
    subplot(nrow,ele,i);
    imagesc(Mt(:,:,i));
    set(gca,'xtick',[],'ytick',[]);
    if(i==1);xlabel('\tau');ylabel('original: \theta');end
    subplot(nrow,ele,i+ele);
    imagesc(prj(:,:,i));
    set(gca,'xtick',[],'ytick',[]);

    if(i==1);xlabel('\tau');ylabel('aligned: \theta');end
    subplot(nrow,ele,ele*2+i);
    imagesc(xc(:,:,i));
    if(i==1);ylabel('reconstruction');end
end
% save(['result/',sample,'/Recon_shift_sc_',num2str(numThetan),'_',num2str(NumElement),'_',num2str(nTau+1),'poiss',num2str(poiss_ratio),'.mat'],'xtot','prj');

