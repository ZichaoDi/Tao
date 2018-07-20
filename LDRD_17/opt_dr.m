%%==Reconstruction of 2D sample with implicit CoR correction
global n_delta N_delta W0 sinoS
global nit maxiter fiter x_iter
N_delta=numThetan;%n_delta/2;
%%==================================================
deltaStar=(cos(theta)-1).*delta0_bar(1:2:n_delta)+sin(theta).*delta0_bar(2:2:n_delta);
x_res=[];
nit_res=[];
aligned=[];
W0=[deltaStar;W(:)];
maxiter=1000;
NF = [0*N; 0*N; 0*N];
Lmap=[];
if(synthetic==0)
    Lmap=sparse(L*Q);
else
    if(ndims(DisR)>=3)
        Mt=squeeze(-log(DisR(:,:,2)./I0'));
        Mt0=Mt;
        Lmap=sparse(squeeze(L_cr(:,:,1)));
        sinoS=squeeze(-log(DisR(:,:,1)./I0'));
    else
        Mt=-log(DisR./I0');
        Mt=Mt-min(Mt(:));
        Lmap=L;
    end
    Lmap=Lmap*Q;
end

for res_step=1:initial_direction
    noise=0*rand(size(Mt0))*res_step*5e-3;
    Mt=Mt0+noise;
    delta=(cos(theta)-1).*x0_opt(1:2:n_delta,res_step)+sin(theta).*x0_opt(2:2:n_delta,res_step);
    x0=[delta;x0_opt(n_delta+1:end,res_step)];
    err0=norm(W0-x0);
    fctn_COR=@(x)sfun_cor_dr(x,full(sinoS'),sparse(Lmap));% on attenuation coefficients miu;
    % fctn_COR=@(x)sfun_pad(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
    low=[floor(-nTau/2*ones(N_delta,1));zeros(N^2*NumElement,1)];
    up=[floor(nTau/2*ones(N_delta,1));inf*ones(N^2*NumElement,1)];
    % low=[-inf*ones(N_delta,1);zeros(prod(m)*NumElement,1)];
    % up=[inf*ones(N_delta,1);inf*ones(prod(m)*NumElement,1)];
    bounds=1;
    if(bounds)
        [xCOR,f,g,ierror] = tn (x0,fctn_COR);%,low,up); % algo='TNbc';
        f1_dr=fiter;
        x1_dr=x_iter;
        xCOR(N_delta+1:end)=Q*xCOR(N_delta+1:end);
        errW_dr(res_step)=norm(W0(N_delta+1:end)-xCOR(N_delta+1:end),'inf');
        f_global(res_step)=f;
        [~,~,alignedSignal]=feval(fctn_COR,xCOR);
        aligned=[aligned,alignedSignal(:)];

        %% ============================================
        % load align86_5;
        % [~,~,xtm1]=feval(fctn_COR,xCOR);
        % % xtm1=max(0,xtm1-0.1);
        % % for n=1:numThetan
        % %     xtm1(n,:)=map1D(xtm1(n,:),[0 1]);
        % % end
        % fctn=@(x)sfun_radon(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
        % W0=W(:);
        % [x,f,g,ierror] = tnbc (x0(N_delta+1:end),fctn,low(N_delta+1:end),up(N_delta+1:end)); % algo='TNbc';
        % W0=[deltaStar;W(:)];
        % %% ============================================

        % err_sino(res_step)=norm(alignedSignal'-sinoS);
        % options = optimoptions('intlinprog','Display','iter');%,'Algorithm','interior-point');
        % [x, f] = intlinprog(fctn_COR,[1:N_delta],[],[],[],[],low,up,options); %algo='fmincon';

        %%============================================
    else
        [xCOR,f,g,ierror] = tn (x0,fctn_COR);
    end
    x_res=[x_res,xCOR];
    nit_res=[nit_res,nit];
end
nrow=4;
figure, 
for i=1:initial_direction, 
    subplot(nrow,initial_direction,i);imagesc(reshape(x_res(N_delta+1:end,i),N,N));%,[min(W(:)),max(W(:))]);
set(gca,'xtick',[],'ytick',[]);
    if(i==1);ylabel({'reconstructed', 'sample'});end
end

subplot(nrow,initial_direction,[1*initial_direction+1:initial_direction*2]);plot(f_global,'ro-'); 
ylabel({'objective', 'value'})
for i=1:initial_direction, 
    subplot(nrow,initial_direction,2*initial_direction+i);imagesc(reshape(aligned(:,i),numThetan,nTau+1));
set(gca,'xtick',[],'ytick',[]);
    if(i==1);ylabel({'corresponding','sinogram','\theta'});xlabel('\tau');end
end
subplot(nrow,initial_direction,[3*initial_direction+1:initial_direction*4]);
plot(errW,'b*-'); ylabel({'reconstruction', 'error'})
