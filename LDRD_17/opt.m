global maxiter W0
rng('default');
x0 = 0*10^(0)*rand(m(1)*m(2)*NumElement,1);
delta=0*delta0';
W0=W(:);
x0=[delta;x0];
W0=[delta0';W0];
err0=norm(W0-x0);
NF = [0*N; 0*N; 0*N];
maxiter=100;
Mt=-log(DisR(:,:,2)./I0);
% for k_step=0
%     delta=delta-k_step*Ddelta.*abs(ceil(eps2))*repmat(dz(1),[numThetan,1]);
%     XTM=aligned(delta,Mt,DetKnot0(1,:),SourceKnot0(1,:),thetan,nTau,numThetan,dTau);
    fctn=@(x)sfun_radonCOR(x,Mt',L_cr(:,:,1));% on attenuation coefficients miu;
    low=[-inf; -inf; zeros(prod(m)*NumElement,1)];
    up=inf*ones(prod(m)*NumElement+2,1);
    bounds=1;
    if(bounds)
        [x,f,g,ierror] = tnbc (x0,fctn,low,up); % algo='TNbc';
        % options = optimoptions('fmincon','Display','iter','GradObj','off','MaxIter',maxiter,'Algorithm','interior-point');% ,'Algorithm','Trust-Region-Reflective'
        % [x, f] = fmincon(fctn,x0,[],[],[],[],low,up,[],options); algo='fmincon';
    else
        [x,f,g,ierror] = tn (x0,fctn);
    end

% end

