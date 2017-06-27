global N_level n_delta maxiter W0
rng('default');
Det=norm(DetKnot0(1,:)-SourceKnot0(1,:))*dTau;
delta=repmat(cr(2,:),[numThetan,1])+pert;
delta0_bar=[];
if(n_delta==2)
    %%======== fixed COR for every projection
    delta0_bar=[((DetKnot0(1,1)-SourceKnot0(1,1))*delta0(2)-(DetKnot0(1,2)-SourceKnot0(1,2))*delta0(1))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*delta0(1)+(DetKnot0(1,2)-SourceKnot0(1,2))*delta0(2))/Det];
elseif(n_delta==2*numThetan)
    %%======== Different COR for each projeciton
    delta0_bar=reshape([((DetKnot0(1,1)-SourceKnot0(1,1))*delta(:,2)-(DetKnot0(1,2)-SourceKnot0(1,2))*delta(:,1))/Det ((DetKnot0(1,1)-SourceKnot0(1,1))*delta(:,1)+(DetKnot0(1,2)-SourceKnot0(1,2))*delta(:,2))/Det]',n_delta,1);
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
res=7;
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,2:4);
initial_direction=size(d0,2);
N_level=[N (N+1)/2];
maxiter=50;
x_res_fine=[];
x_res_coarse=[];
x0_COR=[];
err_sino=[];
f_global=[];
Lmap=[];
Mt0=squeeze(-log(DisR(:,:,1)./I0'));
NF = [0*N_level; 0*N_level; 0*N_level];
Mt=squeeze(-log(DisR(:,:,2)./I0));
% Mt=-log(xrt_shift./I0);

low_CORx=repmat(omega(1),1,n_delta/2);
low_CORy=repmat(omega(3),1,n_delta/2);
low_COR=reshape([low_CORx;low_CORy],n_delta,1)/dTau;
up_CORx=repmat(omega(2),1,n_delta/2);
up_CORy=repmat(omega(4),1,n_delta/2);
up_COR=reshape([up_CORx;up_CORy],n_delta,1)/dTau;
for res_step=1:initial_direction
    Lmap=squeeze(L_cr(:,:,1));
    if(n_delta==3)
        delta=[d0(:,res_step);-res]+1*deltaStar;
    else
        delta=repmat(d0(:,res_step),n_delta/2,1)+1*deltaStar;
    end
    for level=length(N_level):-1:2 
        current_n=N_level(level);
        W0=[deltaStar;downdate(W(:),1)];
        Lmap_coarse=downdate_radon(Lmap,numThetan,nTau);
        Q=diag(1./sum(Lmap_coarse,1));
        Lmap_coarse=Lmap_coarse*Q;
        
        x0=[delta;0*10^(0)*rand(current_n^2*NumElement,1)];
        err0=norm(W0-x0);
        fctn_COR=@(x)sfun_COR(x,full(Mt'),sparse(Lmap_coarse));% on attenuation coefficients miu;

        low=[low_COR;zeros(current_n^2*NumElement,1)];
        up=[up_COR;inf*ones(current_n^2*NumElement,1)];
        bounds=1;
        if(bounds)
            [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
        else
            [x,f,g,ierror] = tn (x0,fctn);
        end

        x_res_coarse=[x_res_coarse,[xCOR(1:n_delta);Q*xCOR(n_delta+1:end)]];
        x0_COR=[x0_COR,x0(1:n_delta)];
        
        x0=[xCOR(1:n_delta);update(1*xCOR(n_delta+1:end),1)];

        current_n=N_level(level-1);

        W0=[deltaStar;W(:)];
        Q=diag(1./sum(Lmap,1));
        Lmap=Lmap*Q;
        err0=norm(W0-x0);
        fctn_COR=@(x)sfun_COR(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
        low=[low_COR;zeros(current_n^2*NumElement,1)];
        up=[up_COR;inf*ones(current_n^2*NumElement,1)];
        bounds=1;
        if(bounds)
            [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
            xCOR(n_delta+1:end)=Q*xCOR(n_delta+1:end);
        else
            [x,f,g,ierror] = tn (x0,fctn);
        end
        [~,~,alignedSignal]=feval(fctn_COR,xCOR);

        % xCOR=x0;
        % [f,g,alignedSignal]=feval(fctn_COR,xCOR);

        f_global(res_step)=f;
        err_sino(res_step)=norm(alignedSignal'-Mt0);
        x_res_fine=[x_res_fine,xCOR];
    end
end

err_fine=sqrt(sum((x_res_fine(n_delta+1:end,:)-repmat(W(:),1,initial_direction)).^2,1));

% figure,
% for i=1:initial_direction, 
%     subplot(initial_direction,2,1+2*(i-1));
%     imagesc(reshape(x_res_coarse(n_delta+1:end,i),N_level(2),N_level(2)));
%     axis xy image; colorbar;
%     subplot(initial_direction,2,2+2*(i-1));
%     imagesc(reshape(x_res_fine(n_delta+1:end,i),N_level(1),N_level(1)));
%     axis xy image; colorbar;
% end
figure, 
nrow=6;
for i=1:initial_direction, 
    subplot(nrow,initial_direction,i);imagesc(reshape(x_res_coarse(n_delta+1:end,i),N_level(2),N_level(2))); colorbar;if(i==1),ylabel('coarse solution');end
    subplot(nrow,initial_direction,i+initial_direction);imagesc(reshape(x_res_fine(n_delta+1:end,i),N,N)); colorbar;if(i==1),ylabel('fine solution');end
    subplot(nrow,initial_direction,i+2*initial_direction);
    plot(x_res_coarse(1:2:n_delta-1,i),x_res_coarse(2:2:n_delta,i),'b.',deltaStar(1:2:n_delta-1),deltaStar(2:2:n_delta),'r*',x0_COR(1:2:end-1,i),x0_COR(2:2:end,i),'g*'); if(i==1),ylabel('CORs-coarse'); legend('reconstructed CORs','true CORs','initial CORs'); end
    for j=1:2:n_delta
        text(x_res_coarse(j,i),x_res_coarse(j+1,i),num2str(j));end

    subplot(nrow,initial_direction,i+3*initial_direction);
    plot(x_res_fine(1:2:n_delta-1,i),x_res_fine(2:2:n_delta,i),'b.',deltaStar(1:2:n_delta-1),deltaStar(2:2:n_delta),'r*');if(i==1),ylabel('CORs-fine');end
end

% subplot(nrow,initial_direction,[1+2*initial_direction:initial_direction*3]);
% h(1)=plot(deltaStar(1:2:n_delta-1),deltaStar(2:2:n_delta),'r*');hold on;
% cmap=hsv(initial_direction);
% for i=1:initial_direction, 
%     h(i+1)=plot(x_res_coarse(1:2:n_delta-1,i),x_res_coarse(2:2:n_delta,i),'.','Color',cmap(i,:));
% end
% hold off; legend(h,'true','first','second');ylabel('Corse Recovered CORs');

% subplot(nrow,initial_direction,[1+3*initial_direction:initial_direction*4]);
% h(1)=plot(deltaStar(1:2:n_delta-1),deltaStar(2:2:n_delta),'r*');hold on;
% cmap=hsv(initial_direction);
% for i=1:initial_direction, 
%     h(i+1)=plot(x_res_fine(1:2:n_delta-1,i),x_res_fine(2:2:n_delta,i),'.','Color',cmap(i,:));
% end
% hold off; legend(h,'true','first','second');ylabel('Fine Recovered CORs');

subplot(nrow,initial_direction,[4*initial_direction+1:initial_direction*5]);
plot(err_sino,'b*-'); ylabel('sinogram error')
subplot(nrow,initial_direction,[5*initial_direction+1:initial_direction*6]);semilogy(f_global,'ro-'); 
ylabel('objective value')


