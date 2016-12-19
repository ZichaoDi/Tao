global n_delta maxiter W0
rng('default');
n_delta=3;
Det=norm(DetKnot0(1,:)-SourceKnot0(1,:))*dTau;
delta0_bar=[((DetKnot0(1,1)-SourceKnot0(1,1))*delta0(2)-(DetKnot0(1,2)-SourceKnot0(1,2))*delta0(1))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*delta0(1)-(DetKnot0(1,2)-SourceKnot0(1,2))*delta0(2))/Det];
if(n_delta==2)
    deltaStar=delta0_bar';%delta0'/dTau;
else
    deltaStar=[delta0_bar';1/2*delta_d0/dTau];
end
res=1;
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
% d0=d0(:,2);
W0=[deltaStar;W(:)];
errW=zeros(size(d0,2),1);
maxiter=50;
for res_step=1:size(d0,2)
    if(n_delta==2)
        delta=0*[d0(:,res_step)]+1*deltaStar;
    else
        delta=[d0(:,res_step);0*deltaStar(3)]+1*deltaStar;
    end
    x0=[delta;0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
    err0=norm(W0-x0);
    NF = [0*N; 0*N; 0*N];
    if(ndims(DisR)>=3)
        Mt=squeeze(-log(DisR(:,:,2)./I0));
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
    fctn=@(x)sfun_radonCOR_sim(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
    low=[-inf*ones(n_delta,1);zeros(prod(m)*NumElement,1)];
    up=inf*ones(prod(m)*NumElement+n_delta,1);
    bounds=1;
    if(bounds)
        [x,f,g,ierror] = tnbc (x0,fctn,low,up); % algo='TNbc';
        % options = optimoptions('fmincon','Display','iter','GradObj','on','MaxIter',maxiter,'Algorithm','sqp');%,'interior-point');% ,'Algorithm'
        % [x, f] = fmincon(fctn,x0,[],[],[],[],low,up,[],options); %algo='fmincon';
    else
        [x,f,g,ierror] = tn (x0,fctn);
    end
    figure(1),
    subplot(size(d0,2),1,res_step);
    imagesc(reshape(x(n_delta+1:end),N,N));
    hold on;
    x_res=[x_res,x];
    errW(res_step)=norm(W0-x);
end
% end

% p=linspace(-20,20,100);
% f=[];
% for i=1:length(p)
%     [f(i)]=feval(fctn,[p(i);p(i);p(i);W(:)]);
% end
