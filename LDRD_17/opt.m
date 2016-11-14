global maxiter W0
rng('default');
x0 = 0*10^(0)*rand(m(1),m(2),NumElement);
x0=x0(:);
delta=[0;0];
x0=[delta;x0];
W0=W(:);
err0=norm(W0(:)-x0(3:end));
NF = [0*N; 0*N; 0*N];
maxiter=100;
Mt=-log(DisR(:,:,2)./I0);

fctn=@(x)sfun_radonCOR(x,Mt,L_cr(:,:,1));% on attenuation coefficients miu;
low=[-inf;-inf;zeros(prod(m)*NumElement,1)];
up=inf*ones(size(x0));
low=low(3:end);
up=up(3:end);
bounds=1;
if(bounds)
    [x,f,g,ierror] = tnbc (x0(3:end),fctn,low,up); % algo='TNbc';
    % options = optimoptions('fmincon','Display','iter','GradObj','off','MaxIter',maxiter,'Algorithm','interior-point');% ,'Algorithm','Trust-Region-Reflective'
    % [x, f] = fmincon(fctn,x0,[],[],[],[],low,up,[],options); algo='fmincon';
else
    [x,f,g,ierror] = tn (x0,fctn);
end



