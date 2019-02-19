%% Reconstruction of 2D sample with explicit CoR correction
global n_delta maxiter W0 sinoS
global x_iter fiter nit
global reg
rng('default');
Det=norm(DetKnot0(1,:)-SourceKnot0(1,:))*dTau;
delta0_bar=[];
if(n_delta==2)
    %%======== fixed COR for every projection
    delta0_bar=[((DetKnot0(1,1)-SourceKnot0(1,1))*delta0(1,2)-(DetKnot0(1,2)-SourceKnot0(1,2))*delta0(1,1))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*delta0(1,1)+(DetKnot0(1,2)-SourceKnot0(1,2))*delta0(1,2))/Det]';
elseif(n_delta==2*numThetan)
    %%======== Different COR for each projeciton
    delta0_bar=reshape([((DetKnot0(1,1)-SourceKnot0(1,1))*delta0(:,2)-(DetKnot0(1,2)-SourceKnot0(1,2))*delta0(:,1))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*delta0(:,1)+(DetKnot0(1,2)-SourceKnot0(1,2))*delta0(:,2))/Det]',n_delta,1);
elseif(n_delta==4)
    %%======== single COR for half of the projections
    delta0_bar(1:2)=[((DetKnot0(1,1)-SourceKnot0(1,1))*(cr(2,2)+pert1)-(DetKnot0(1,2)-SourceKnot0(1,2))*(cr(2,1)+pert1))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*(cr(2,1)+pert1)-(DetKnot0(1,2)-SourceKnot0(1,2))*(cr(2,2)+pert1))/Det];
    delta0_bar(3:4)=[((DetKnot0(1,1)-SourceKnot0(1,1))*(cr(2,2)+pert2)-(DetKnot0(1,2)-SourceKnot0(1,2))*(cr(2,1)+pert2))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*(cr(2,1)+pert2)-(DetKnot0(1,2)-SourceKnot0(1,2))*(cr(2,2)+pert2))/Det];
    delta0_bar=delta0_bar';
end
%%==================================================
if(n_delta==3)
    deltaStar=[delta0_bar';1/2*delta_d0/dTau];
else
    deltaStar=delta0_bar;%delta0'/dTau;
end
Lmap=[];
if(synthetic==0)
    deltaStar=0*deltaStar;
    NumElement=1;
    ele_ind=4;
    % Z=Z(ele_ind);
    % W=W(:,:,ele_ind);
    Mt=XRF_raw_tot(:,:,ele_ind,86);
    Mt=Mt./max(Mt(:));%%==normalize data;
    Mt=Mt-min(Mt(:));%%==normalize data;
    figure, imagesc(Mt)
    Mt0=Mt;
    % for i=1:numThetan,
    %     % Mt(:,i)=medfilt1(Mt(:,i),3);
    %     Mt(:,i)=max(0,Mt(:,i)-0.05);%min(Mt(:,i));
    % end
    sinoS=Mt;
    % Q=sparse(diag(1./sum(L,1)));
    Q=speye(size(L,2));
    Lmap=sparse(L*Q);
else
    if(ndims(DisR)>=3)
        Mt=squeeze(-log(DisR(:,:,2)./I0'));%/scale;
        Mt0=Mt;
        Lmap=sparse(squeeze(L_cr(:,:,1)));
        sinoS=squeeze(-log(DisR(:,:,1)./I0'));
    else
        Mt=-log(DisR./I0');
        Lmap=sparse(L);
    end
    Q=sparse(diag(1./sum(Lmap,1).^(1/2)));
    Q=1e0*speye(size(Q));
    Lmap=Lmap*Q;
end
res=7;
x_res=[];
nit_res=[];
aligned=[];
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,[1]);% 3 5 7 9]);
% d0=repmat(d0(:,1),1,10);
initial_direction=size(d0,2);
W0=[deltaStar;W(:)];
errW=zeros(size(d0,2),1);
f_global=zeros(size(d0,2),1);
x0_opt=[];
maxiter=100;
NF = [0*N; 0*N; 0*N];
% noise_level=[0.04 0.09 0.13 0.18 0.22];
% Mt_noise=zeros(numThetan,nTau+1,length(noise_level));
% rng(1);
% for i=1:5
%     noise=rand(size(Mt))*res_step*noise_level(i);
%     Mt_noise(:,:,i)=Mt+noise;
% end
% save('noise_1cor.mat','Mt_noise');
% return;

%%====bound shift parameter by the boundary of the sinogram
low_CORx=repmat(omega(1),1,n_delta/2);
low_CORy=repmat(omega(3),1,n_delta/2);
low_COR=reshape([low_CORx;low_CORy],n_delta,1)/dTau;
up_CORx=repmat(omega(2),1,n_delta/2);
up_CORy=repmat(omega(4),1,n_delta/2);
up_COR=reshape([up_CORx;up_CORy],n_delta,1)/dTau;
% low=[low_COR;zeros(prod(m)*NumElement,1)];
% up=[up_COR;inf*ones(prod(m)*NumElement,1)];
%%====================================================
low=[-nTau/2.*ones(size(low_COR));zeros(prod(m)*NumElement,1)];
up=[nTau/2.*ones(size(up_COR));inf*ones(prod(m)*NumElement,1)];
% %% ============================================
bounds=1;
if(synthetic==0)
    fctn=@(x)sfun_radon(x,Mt,sparse(Lmap));% on attenuation coefficients miu;
else
    fctn=@(x)sfun_radon(x,Mt,squeeze(L_cr(:,:,2)));% on attenuation coefficients miu;
end;
W0=W(:);
err0=norm(W0);
[x_base,f,g,ierror] = tnbc(ones(N^2,1),fctn,low(n_delta+1:end),up(n_delta+1:end)); % algo='TNbc';
% mt=0.*Mt; 
% mirror_shift=16;
% mt(:,mirror_shift:end)=Mt(:,1:end-mirror_shift+1);
% [x_mirror,f,g,ierror] = tnbc(ones(N^2,1),fctn,low(n_delta+1:end),up(n_delta+1:end)); % algo='TNbc';
% save('fig7_extra.mat','x_base','x_mirror');
f0=fiter;xiter0=x_iter;
W0=[deltaStar;W(:)];

% %% ============================================
    
for res_step=1:initial_direction
    noise=0;%rand(size(Mt))*res_step*noise_level(i);
    Mt=Mt+noise;
    per=0;
    if(n_delta==3)
        delta=[d0(:,res_step);-res]+0*deltaStar;
    else
        delta=repmat(d0(:,res_step),n_delta/2,1)+per*deltaStar;
    end
    x0=[delta;0*rand(m(1)*m(2)*NumElement,1)];
    x0_opt(:,res_step)=x0;
    err0=norm(W0-x0);
    fctn_COR=@(x)sfun_cor(x,full(Mt),sparse(Lmap));% on attenuation coefficients miu;
    if(bounds)
        [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
        xCOR(n_delta+1:end)=Q*xCOR(n_delta+1:end);
        errW(res_step)=norm(W0(n_delta+1:end)-xCOR(n_delta+1:end));
        f_global(res_step)=f;
        % load aligned86_5_paunesku;
        [~,~,alignedSignal]=feval(fctn_COR,xCOR);
        err_sino(res_step)=norm(alignedSignal-sinoS);
    else
        [xCOR,f,g,ierror] = tn (x0,fctn_COR);
    end
    f1=fiter;xiter1=x_iter;
    x_res=[x_res,xCOR];
    nit_res=[nit_res,nit];
    aligned=[aligned,alignedSignal(:)];
end
nrow=4;
figure,
for i=1:initial_direction,
    subplot(nrow,initial_direction,i);imagesc(reshape(x_res(n_delta+1:end,i),N,N));%,[min(W(:)),max(W(:))]);
    set(gca,'xtick',[],'ytick',[]);
    if(i==1);ylabel({'reconstructed', 'sample'});end
end

% subplot(nrow,initial_direction,[1+initial_direction:initial_direction*2]);
% plot(deltaStar(1:2:n_delta),deltaStar(2:2:n_delta),'r*');hold on;
% legends{1}='true';
% cmap=hsv(initial_direction+1);
% plot(x0(1:2:n_delta),x0(2:2:n_delta),'o','Color',cmap(1,:));
% legends{2}='initial';
% for i=1:initial_direction,
%     plot(x_res(1:2:n_delta-1,i),x_res(2:2:n_delta,i),'.','Color',cmap(i+1,:));
% end
% legends{3}='recovered';
% hold off; legend(legends);ylabel('Recovered CORs');
subplot(nrow,initial_direction,[1*initial_direction+1:initial_direction*2]);plot(f_global,'ro-');
ylabel({'objective', 'value'})
for i=1:initial_direction,
    subplot(nrow,initial_direction,2*initial_direction+i);imagesc(reshape(aligned(:,i),numThetan,nTau+1));
    set(gca,'xtick',[],'ytick',[]);
    if(i==1);ylabel({'corresponding','sinogram','\theta'});xlabel('\tau');end
end
subplot(nrow,initial_direction,[3*initial_direction+1:initial_direction*4]);
plot(errW,'b*-'); ylabel({'reconstruction', 'error'})

