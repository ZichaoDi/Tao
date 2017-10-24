global n_delta N_delta maxiter W0 sinoS
N_delta=numThetan;%n_delta/2;
%%==================================================
deltaStar=(cos(theta)-1).*delta0_bar(1:2:n_delta)+sin(theta).*delta0_bar(2:2:n_delta);
% res=2;
x_res=[];
aligned=[];
W0=[deltaStar;W(:)];
maxiter=150;
NF = [0*N; 0*N; 0*N];
Lmap=[];
if(synthetic==0)
    Lmap=sparse(L*Q);
else
    if(ndims(DisR)>=3)
        Mt=squeeze(-log(DisR(:,:,2)./I0'));
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
    delta=(cos(theta)-1).*x0_opt(1:2:n_delta,res_step)+sin(theta).*x0_opt(2:2:n_delta,res_step);
    x0=[delta;x0_opt(n_delta+1:end,res_step)];
    err0=norm(W0-x0);
    fctn_COR=@(x)sfun_cor_dr(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
    low=[floor(-nTau/2*ones(N_delta,1));zeros(prod(m)*NumElement,1)];
    up=[floor(nTau/2*ones(N_delta,1));inf*ones(prod(m)*NumElement,1)];
    % low=[-inf*ones(N_delta,1);zeros(prod(m)*NumElement,1)];
    % up=[inf*ones(N_delta,1);inf*ones(prod(m)*NumElement,1)];
    bounds=1;
    if(bounds)
        [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
        xCOR(N_delta+1:end)=Q*xCOR(N_delta+1:end);
        errW_dr(res_step)=norm(W0(N_delta+1:end)-xCOR(N_delta+1:end));
        f_global(res_step)=f;
        [~,~,alignedSignal]=feval(fctn_COR,xCOR);
        aligned=[aligned,alignedSignal(:)];

        %% ============================================
        % load align86_5;
        [~,~,xtm1]=feval(fctn_COR,xCOR);
        % xtm1=max(0,xtm1-0.1);
        % for n=1:numThetan
        %     xtm1(n,:)=map1D(xtm1(n,:),[0 1]);
        % end
        fctn=@(x)sfun_radon(x,full(xtm1),sparse(Lmap));% on attenuation coefficients miu;
        W0=W(:);
        [x,f,g,ierror] = tnbc (x0(N_delta+1:end),fctn,low(N_delta+1:end),up(N_delta+1:end)); % algo='TNbc';
        W0=[deltaStar;W(:)];
        % %% ============================================

        % err_sino(res_step)=norm(alignedSignal'-sinoS);
        % options = optimoptions('intlinprog','Display','iter');%,'Algorithm','interior-point');
        % [x, f] = intlinprog(fctn_COR,[1:N_delta],[],[],[],[],low,up,options); %algo='fmincon';

        %%============================================
    else
        [xCOR,f,g,ierror] = tn (x0,fctn_COR);
    end
    x_res=[x_res,xCOR];
end
figure,
nrow=4;
for i=1:initial_direction, 
    subplot(nrow,initial_direction,i);imagesc(reshape(x_res(N_delta+1:end,i),N,N));
    if(i==1);ylabel('reconstructed sample');end
end
% subplot(nrow,initial_direction,[1+initial_direction:initial_direction*2]);
% plot(deltaStar,zeros(N_delta,1),'r*');hold on;
% legends{1}='true';
% cmap=hsv(initial_direction+1);
% plot(x0(1:N_delta),zeros(N_delta,1),'o','Color',cmap(1,:));
% legends{2}='initial';
% for i=1:initial_direction, 
%     plot(x_res(1:N_delta,i),zeros(N_delta,1),'.','Color',cmap(i+1,:));
% end
% legends{3}='recovered';
% hold off; legend(legends);ylabel('Recovered CORs');

subplot(nrow,initial_direction,[1*initial_direction+1:initial_direction*2]);
plot(errW_dr,'b*-'); ylabel('reconstruction error')
subplot(nrow,initial_direction,[2*initial_direction+1:initial_direction*3]);plot(f_global,'ro-'); 
ylabel('objective value')
for i=1:initial_direction, 
    subplot(nrow,initial_direction,3*initial_direction+i);imagesc(reshape(aligned(:,i),numThetan,nTau+1));
end

