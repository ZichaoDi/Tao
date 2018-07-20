global n_delta N_delta maxiter W0 sinoS
N_delta=numThetan;%n_delta/2;
%%==================================================
ntot=N_delta+N^2;
deltaStar=zeros(N_delta,1);
W0=[deltaStar;W(:)];
maxiter=200;
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
end

x0=zeros(ntot,1);
err0=norm(W0-x0);
% fctn_COR=@(x)sfun_cor_dr(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
fctn_COR=@(x)sfun_pad(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
low=[floor(-nTau/2*ones(N_delta,1));zeros(N^2*NumElement,1)];
up=[floor(nTau/2*ones(N_delta,1));inf*ones(N^2*NumElement,1)];
% low=[-inf*ones(N_delta,1);zeros(prod(m)*NumElement,1)];
% up=[inf*ones(N_delta,1);inf*ones(prod(m)*NumElement,1)];
bounds=1;
[xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
[~,~,alignedSignal]=feval(fctn_COR,xCOR);
%% ============================================
fctn=@(x)sfun_radon(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
W0=W(:);
[x,f,g,ierror] = tnbc (x0(N_delta+1:end),fctn,low(N_delta+1:end),up(N_delta+1:end)); % algo='TNbc';
W0=[deltaStar;W(:)];
% %% ============================================
