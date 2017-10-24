global N_delta maxiter W0 sinoS
global xiter fiter ErrIter
global alpha

ReconAttenu = 0; % 0: Recover W; 1: Recover miu
slice=1;
do_setup;
N_delta=numThetan;
rng('default');
res=3;
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,3);
initial_direction=size(d0,2);
deltaStar=zeros(N_delta,1);
W0=[deltaStar;W(:)];
maxiter=200;
recon=sparse(N^2*NumElement+N_delta,nslice);
for slice=1:nslice
    setup_paunesku;
    NF = [0*N; 0*N; 0*N];
    Mt=XRF_decom(:,:,:);
    for ele=1:NumElement
        a=Mt(:,:,ele);alpha(ele)=std(a(:));
    end
    alpha=alpha./sum(alpha);
    Mt=Mt./max(Mt(:))+abs(min(Mt(:)));%%==normalize data;
    sinoS=reshape(Mt,numThetan*(nTau+1),NumElement);
    Lmap=sparse(L);
    for res_step=1:initial_direction
        delta=repmat(d0(:,res_step),n_delta/2,1);
        x0 = [(cos(theta)-1)*delta(1)+sin(theta)*delta(2);0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
        err0=norm(W0-x0);
        fctn_COR=@(x)sfun_shift_mc(x,Mt,Lmap);% on attenuation coefficients miu;
        low=[-nTau/2*ones(N_delta,1);zeros(prod(m)*NumElement,1)];
        up=[nTau/2*ones(N_delta,1);inf*ones(prod(m)*NumElement,1)];
        bounds=1;
        if(bounds)
            [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
            errW(res_step)=norm(W0(n_delta+1:end)-xCOR(n_delta+1:end));
            f_global(res_step)=f;
            [~,~,alignedSignal]=feval(fctn_COR,xCOR);
            %% ============================================
            % x1=zeros(N^2,NumElement);
            % for ele=1:NumElement
            %     xtm1=reshape(alignedSignal,numThetan,nTau+1,NumElement);
            %     fctn=@(x)sfun_radon(x,full(xtm1(:,:,ele)),sparse(Lmap));% on attenuation coefficients miu;
            %     W0=zeros(N^2,1);
            %     [x,f,g,ierror] = tnbc (W0,fctn,zeros(N^2,1),inf*ones(N^2,1)); % algo='TNbc';
            %     x1(:,ele)=x;
            % end
            %% ============================================
        else
            [xCOR,f,g,ierror] = tn (x0,fctn_COR);
        end
    end
    recon(:,slice)=xCOR;
end
save(['result/miller/Recon_shift_mc_weighted.mat'],'recon');

% nrow=3;
% figure,
% xc=reshape(x_res(N_delta+1:end,:),[N,N,NumElement,initial_direction]);
% for ele=1:NumElement
%     subplot(nrow,NumElement,ele);
%     imagesc(XRF_decom(:,:,ele));
%     if(ele==1);ylabel('original');end
%     subplot(nrow,NumElement,NumElement+ele);
%     align=reshape(alignedSignal,numThetan*(nTau+1),NumElement);
%     imagesc(reshape(align(:,ele),numThetan,nTau+1));
%     if(ele==1);ylabel('aligned');end
%     subplot(nrow,NumElement,NumElement*2+ele);
%     imagesc(xc(:,:,ele));
%     if(ele==1);ylabel('reconstruction');end
% end
