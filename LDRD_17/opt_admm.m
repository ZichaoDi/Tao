global n_delta maxiter W0 sinoS
global xiter fiter ErrIter N_delta
rng('default');
Det=norm(DetKnot0(1,:)-SourceKnot0(1,:))*dTau;
deltaStar=[];
if(n_delta==2*numThetan)
    %%======== Different COR for each projeciton
    deltaStar=reshape([((DetKnot0(1,1)-SourceKnot0(1,1))*delta0(:,2)-(DetKnot0(1,2)-SourceKnot0(1,2))*delta0(:,1))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*delta0(:,1)+(DetKnot0(1,2)-SourceKnot0(1,2))*delta0(:,2))/Det]',n_delta,1);
end
delta0_bar=deltaStar;
%%==================================================
res=2;
x_res=[];
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,6);
initial_direction=size(d0,2);
W0=[deltaStar;W(:)];
errW=[];
f_global=zeros(size(d0,2),1);
maxiter=50;
NF = [0*N; 0*N; 0*N];
Lmap=[];
if(synthetic==0)
    Mt=XRF_decom(:,:,1)';
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
    Q=sparse(diag(1./sum(Lmap,1).^2));
    Q=speye(size(Q));
    Lmap=Lmap*Q;
end
for res_step=1:initial_direction
    delta=repmat(d0(:,res_step),n_delta/2,1)+1*deltaStar;
    x0=[1*delta;0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
    x0_opt=x0;
    err0=norm(W0-x0);
    low_x=[zeros(prod(m)*NumElement,1)];
    up_x=[inf*ones(prod(m)*NumElement,1)];
    alignedSignal=full(Mt');
    xstar=x0(n_delta+1:end);
    xCOR=x0(1:n_delta);
    % figure,
    for inter_k=1:10
        W0=W(:);
        fctn=@(x)sfun_radon(x,alignedSignal,sparse(Lmap));% on attenuation coefficients miu;
        xold=xstar;
        % [xstar,f,g,ierror] = tn(xstar,fctn); 
        [xstar,f,g,ierror] = tnbc (xstar,fctn,low_x,up_x); 
        fctn_cor=@(x)sfun_cor(x,xstar,alignedSignal,sparse(Lmap));% on attenuation coefficients miu;
        W0=deltaStar;
        gx=g;
        [f,g,xtm]=feval(fctn_cor,xCOR);
        figure, subplot(1,2,1),plot(g);subplot(1,2,2),plot(gx);pause;
        p=-1e84*g;
        [xCOR, f, g, ~,~, alpha] = lin1 (p, xCOR, f, 1, g, fctn_cor);
        alpha
        % [xCOR,f,g,ierror] = tn (xCOR,fctn_cor); % algo='TNbc';
        [~,~,alignedSignal]=feval(fctn_cor,xCOR);
        %%============================================
        errW(inter_k)=norm(W(:)-xstar);
        % plot(inter_k,errW(inter_k),'r.-');hold on; drawnow;
        if(errW(inter_k)<1e-3)
            break;
            disp('converge')
        end
    end
    f_global(res_step)=f;
    err_sino(res_step)=norm(alignedSignal'-sinoS);
    x_res=[x_res,xCOR];
end
nrow=3;
figure, 
for i=1:initial_direction, 
    subplot(nrow,initial_direction,i);imagesc(reshape(xstar,N,N));
    if(i==1);ylabel('reconstructed sample');end
end

subplot(nrow,initial_direction,[1+initial_direction:initial_direction*2]);
plot(deltaStar(1:2:n_delta),deltaStar(2:2:n_delta),'r*');hold on;
legends{1}='true';
cmap=hsv(initial_direction+1);
plot(x0(1:2:n_delta),x0(2:2:n_delta),'o','Color',cmap(1,:));
legends{2}='initial';
for i=1:initial_direction, 
    plot(x_res(1:2:n_delta,i),x_res(2:2:n_delta,i),'.','Color',cmap(i+1,:));
end
legends{3}='recovered';
hold off; legend(legends);ylabel('Recovered CORs');

% subplot(nrow,initial_direction,[2*initial_direction+1:initial_direction*3]);
% plot(err_sino,'b*-'); ylabel('sinogram error')
subplot(nrow,initial_direction,[2*initial_direction+1:initial_direction*3]);semilogy(f_global,'ro-'); 
ylabel('objective value')
% axis_limits=axis; axis([axis_limits-[0 0 abs(axis_limits(3)/6) 0]]);
% end
