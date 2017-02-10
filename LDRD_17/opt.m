global n_delta maxiter W0
rng('default');
Det=norm(DetKnot0(1,:)-SourceKnot0(1,:))*dTau;
delta=repmat(cr(2,:),[numThetan,1])+pert;
if(n_delta==2)
    %%======== fixed COR for every projection
    delta0_bar=[((DetKnot0(1,1)-SourceKnot0(1,1))*delta0(2)-(DetKnot0(1,2)-SourceKnot0(1,2))*delta0(1))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*delta0(1)+(DetKnot0(1,2)-SourceKnot0(1,2))*delta0(2))/Det];
elseif(n_delta==2*numThetan)
    %%======== Different COR for each projeciton
    % for n=1:numThetan
    %     delta0_bar(2*n-1:2*n)=[((DetKnot0(1,1)-SourceKnot0(1,1))*delta(n,2)-(DetKnot0(1,2)-SourceKnot0(1,2))*delta(n,1))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*delta(n,1)+(DetKnot0(1,2)-SourceKnot0(1,2))*delta(n,2))/Det];
    % end
    delta0_bar=reshape([((DetKnot0(1,1)-SourceKnot0(1,1))*delta(:,2)-(DetKnot0(1,2)-SourceKnot0(1,2))*delta(:,1))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*delta(:,1)+(DetKnot0(1,2)-SourceKnot0(1,2))*delta(:,2))/Det]',n_delta,1);
elseif(n_delta==4)
    %%======== single COR for half of the projections
    delta0_bar(1:2)=[((DetKnot0(1,1)-SourceKnot0(1,1))*(cr(2,2)+pert1)-(DetKnot0(1,2)-SourceKnot0(1,2))*(cr(2,1)+pert1))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*(cr(2,1)+pert1)-(DetKnot0(1,2)-SourceKnot0(1,2))*(cr(2,2)+pert1))/Det];
    delta0_bar(3:4)=[((DetKnot0(1,1)-SourceKnot0(1,1))*(cr(2,2)+pert2)-(DetKnot0(1,2)-SourceKnot0(1,2))*(cr(2,1)+pert2))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*(cr(2,1)+pert2)-(DetKnot0(1,2)-SourceKnot0(1,2))*(cr(2,2)+pert2))/Det];
end
%%==================================================
if(n_delta==3)
    deltaStar=[delta0_bar';1/2*delta_d0/dTau];
else
    deltaStar=delta0_bar;%delta0'/dTau;
end
res=0;
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,2:end);
W0=[deltaStar;W(:)];
errW=zeros(size(d0,2),1);
maxiter=50;
for res_step=1:size(d0,2)
    if(n_delta==3)
        delta=[d0(:,res_step);-res]+1*deltaStar;
    else
        delta=res*ones(n_delta,1)+1*deltaStar;
    end
    x0=[delta;0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
    err0=norm(W0-x0);
    NF = [0*N; 0*N; 0*N];
    if(ndims(DisR)>=3)
        Mt=squeeze(-log(DisR(:,:,2)./I0));
        % Mt=-log(xrt_shift./I0);
        Lmap=squeeze(L_cr(:,:,1));
    else
        Mt=-log(DisR./I0');
        Mt=Mt-min(Mt(:));
        % Mt=XRF_decom(:,:,1)';
        Lmap=L;
    end
    % for k_step=0
    %     delta=delta-k_step*Ddelta.*abs(ceil(eps2))*repmat(dz(1),[numThetan,1]);
    %     XTM=aligned(delta,Mt,DetKnot0(1,:),SourceKnot0(1,:),thetan,nTau,numThetan,dTau);
    fctn_COR=@(x)sfun_COR(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
    fctn=@(x)sfun_radon(x,full(squeeze(-log(DisR(:,:,1)./I0))'),sparse(Lmap));% on attenuation coefficients miu;
    low=[-inf*ones(n_delta,1);zeros(prod(m)*NumElement,1)];
    up=inf*ones(prod(m)*NumElement+n_delta,1);
    bounds=1;
    if(bounds)
        [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
        W0=W(:);
        [x,f,g,ierror] = tnbc (x0(n_delta+1:end),fctn,low(n_delta+1:end),up(n_delta+1:end)); % algo='TNbc';
        % options = optimoptions('fmincon','Display','iter','GradObj','on','MaxIter',maxiter,'Algorithm','sqp');%,'interior-point');% ,'Algorithm'
        % [x, f] = fmincon(fctn,x0,[],[],[],[],low,up,[],options); %algo='fmincon';
    else
        [x,f,g,ierror] = tn (x0,fctn);
    end
    figure(4),
    subplot(size(d0,2),1,res_step);
    imagesc(reshape(xCOR(n_delta+1:end),N,N));
    axis xy image;
    hold on;
    x_res=[x_res,xCOR];
    errW(res_step)=norm(W0-xCOR);
end
% end

% p=linspace(-20,20,100);
% f=[];
% for i=1:length(p)
%     [f(i)]=feval(fctn,[p(i);p(i);p(i);W(:)]);
% end
