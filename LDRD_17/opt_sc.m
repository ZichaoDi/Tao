global N_delta maxiter W0 sinoS
global xiter fiter ErrIter
global alpha lambda
global reg_str

ReconAttenu = 0; % 0: Recover W; 1: Recover miu
slice=1;
synthetic=0;
ele_ind=1;
do_setup;
N_delta=numThetan;
rng('default');
res=3;
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,1);
initial_direction=size(d0,2);
deltaStar=zeros(N_delta,1);
W0=[deltaStar;W(:)];
maxiter=100;
aligned3D=zeros(numThetan*(nTau+1),NumElement,nslice);
recon=sparse(N^2*NumElement+N_delta,nslice);
for ele=1:size(XRF_raw_tot,3);
    for slice= 1:nslice
        NF = [0*N; 0*N; 0*N];
        Mt(Mt<0)=0;
        MtS=Mt(:,:,ele);%XRF_raw_tot(:,:,ele,slice);
        for res_step=1:initial_direction
            delta=repmat(d0(:,res_step),n_delta/2,1);
            x0 = [(cos(theta)-1)*delta(1)+sin(theta)*delta(2);0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
            err0=norm(W0-x0);
            alpha=1;
            fctn_COR=@(x)sfun_shift_mc(x,MtS,Lmap);% on attenuation coefficients miu;
            low=[-nTau/2*ones(N_delta,1);zeros(prod(m)*NumElement,1)];
            up=[nTau/2*ones(N_delta,1);inf*ones(prod(m)*NumElement,1)];
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
                    %% ============================================
                else
                    [xCOR,f,g,ierror] = tn (x0,fctn_COR);
                end
            end
        end
        recon(:,slice)=xCOR;
        aligned3D(:,:,slice)=alignedSignal(:);
    end
prj(:,:,ele)=reshape(alignedSignal,numThetan,nTau+1,NumElement);
xtot(:,ele)=xCOR;
nrow=3;
figure,
xc=reshape(xCOR(N_delta+1:end,:),[N,N,NumElement,initial_direction]);
    subplot(nrow,NumElement,1);
    imagesc(XRF_raw_tot(:,:,ele,slice));
    title(Element{Z});
    set(gca,'xtick',[],'ytick',[]);
    if(ele==1);xlabel('\tau');ylabel('original: \theta');end
    subplot(nrow,NumElement,2);
    align=reshape(alignedSignal,numThetan*(nTau+1),NumElement);
    imagesc(reshape(align,numThetan,nTau+1));
    set(gca,'xtick',[],'ytick',[]);
    
    if(ele==1);xlabel('\tau');ylabel('aligned: \theta');end
    subplot(nrow,NumElement,3);
    imagesc(xc(:,:));
    if(ele==1);ylabel('reconstruction');end

end

% save(['result/',sample,'/Recon_shift_sc_noise',num2str(numThetan),'.mat'],'xtot','prj');
