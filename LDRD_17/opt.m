global n_delta maxiter W0 sinoS
global xiter fiter ErrIter
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
res=2;
x_res=[];
aligned=[];
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,2);
initial_direction=size(d0,2);
W0=[deltaStar;W(:)];
errW=zeros(size(d0,2),1);
f_global=zeros(size(d0,2),1);
maxiter=250;
NF = [0*N; 0*N; 0*N];
Lmap=[];
if(synthetic==0)
    deltaStar=0*deltaStar;
    Mt=XRF_decom(:,:,3)';
    Q=sparse(diag(1./sum(L,1)));
    Q=speye(size(Q));
    Lmap=sparse(L*Q);
else
    if(ndims(DisR)>=3)
        Mt=squeeze(-log(DisR(:,:,2)./I0'));
        Lmap=sparse(squeeze(L_cr(:,:,1)));
        sinoS=squeeze(-log(DisR(:,:,1)./I0'));
    else
        Mt=-log(DisR./I0');
        Mt=Mt-min(Mt(:));
        Lmap=sparse(L);
    end
    Q=sparse(diag(1./sum(Lmap,1).^1));
    Q=1e0*speye(size(Q));
    Lmap=Lmap*Q;
end
low_CORx=repmat(omega(1),1,n_delta/2);
low_CORy=repmat(omega(3),1,n_delta/2);
low_COR=reshape([low_CORx;low_CORy],n_delta,1)/dTau;
up_CORx=repmat(omega(2),1,n_delta/2);
up_CORy=repmat(omega(4),1,n_delta/2);
up_COR=reshape([up_CORx;up_CORy],n_delta,1)/dTau;
    % low=[-inf.*ones(size(low_COR));zeros(prod(m)*NumElement,1)];
    % up=[inf.*ones(size(up_COR));inf*ones(prod(m)*NumElement,1)];
    low=[low_COR;zeros(prod(m)*NumElement,1)];
    up=[up_COR;inf*ones(prod(m)*NumElement,1)];
for res_step=1:initial_direction
    if(n_delta==3)
        delta=[d0(:,res_step);-res]+1*deltaStar;
    else
        delta=repmat(d0(:,res_step),n_delta/2,1)+1*deltaStar;
    end
    x0=[delta;0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
    x0_opt=x0;
    err0=norm(W0-x0);
    fctn_COR=@(x)sfun_COR(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
    fctn=@(x)sfun_radon(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
    bounds=1;
    if(bounds)
        [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
        xCOR(n_delta+1:end)=Q*xCOR(n_delta+1:end);
        errW(res_step)=norm(W0(n_delta+1:end)-xCOR(n_delta+1:end));
        f_global(res_step)=f;
        [~,~,alignedSignal]=feval(fctn_COR,xCOR);
        % err_sino(res_step)=norm(alignedSignal'-sinoS);
        %%============================================
        % W0=W(:);
        % [x,f,g,ierror] = tnbc (x0(n_delta+1:end),fctn,low(n_delta+1:end),up(n_delta+1:end)); % algo='TNbc';
        % W0=[deltaStar;W(:)];
        %%============================================
        % options = optimoptions('fmincon','Display','iter','GradObj','on','MaxIter',maxiter);%,'Algorithm','interior-point');
        % [x, f] = fmincon(fctn_COR,x0,[],[],[],[],low,up,[],options); %algo='fmincon';
    else
        [xCOR,f,g,ierror] = tn (x0,fctn_COR);
    end
    x_res=[x_res,xCOR];
    aligned=[aligned,alignedSignal(:)];
end
nrow=4;
figure, 
for i=1:initial_direction, 
    subplot(nrow,initial_direction,i);imagesc(reshape(x_res(n_delta+1:end,i),N,N));
    if(i==1);ylabel('reconstructed sample');end
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
subplot(nrow,initial_direction,[1*initial_direction+1:initial_direction*2]);
plot(errW,'b*-'); ylabel('reconstruction error')
subplot(nrow,initial_direction,[2*initial_direction+1:initial_direction*3]);semilogy(f_global,'ro-'); 
ylabel('objective value')
for i=1:initial_direction, 
    subplot(nrow,initial_direction,3*initial_direction+i);imagesc(reshape(aligned(:,i),numThetan,nTau+1));
end
