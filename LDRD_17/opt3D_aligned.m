global N_delta maxiter W0 sinoS
global xiter fiter ErrIter

slice=1;
do_setup;
N_delta=numThetan;
rng('default');
res=3;
x_res=[];
aligned=[];
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,5);
initial_direction=size(d0,2);
deltaStar=zeros(N_delta,1);
W0=W(:);
maxiter=150;
for element_ind=1:6
    recon=sparse(N^2,131);
    for slice=1:131
        setup_paunesku;
        NF = [0*N; 0*N; 0*N];
        Mt=XRF_decom(:,:,element_ind)';
        Mt=Mt./max(Mt(:));
        sinoS=Mt;
        Lmap=sparse(L);
        low=[-inf.*ones(N_delta,1);zeros(prod(m)*NumElement,1)];
        up=[inf.*ones(N_delta,1);inf*ones(prod(m)*NumElement,1)];
        for res_step=1:initial_direction
            delta=(cos(theta)-1).*d0(1,res_step)+sin(theta).*d0(2,res_step);
            x0= [delta;0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
            err0=norm(x0-x0);
            fctn_COR=@(x)sfun_cor_dr(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
            load align85_4
            x0(1:N_delta)=align85_4;
            [~,~,aligned_xtm]=feval(fctn_COR,x0);
            fctn=@(x)sfun_radon(x,full(aligned_xtm),sparse(Lmap));% on attenuation coefficients miu;
            bounds=1;
            if(bounds)
                [x,f,g,ierror] = tnbc (x0(N_delta+1:end),fctn,low(N_delta+1:end),up(N_delta+1:end)); % algo='TNbc';
            else
                [xCOR,f,g,ierror] = tn (x0,fctn_COR);
            end
            x_res=[x_res,x];
            aligned=[aligned,aligned_xtm(:)];
        end
        recon(:,slice)=x;
    end
    save(['result/paunesku/recon_aligned_reduced_normalized',num2str(element_ind),'.mat'],'recon');
end
% nrow=4;
% figure, 
% for i=1:initial_direction, 
%     subplot(nrow,initial_direction,i);imagesc(reshape(x_res(N_delta+1:end,i),N,N));
%     if(i==1);ylabel('reconstructed sample');end
% end
% 
% subplot(nrow,initial_direction,[1*initial_direction+1:initial_direction*2]);
% plot(errW,'b*-'); ylabel('reconstruction error')
% subplot(nrow,initial_direction,[2*initial_direction+1:initial_direction*3]);plot(f_global,'ro-'); 
% ylabel('objective value')
% for i=1:initial_direction, 
%     subplot(nrow,initial_direction,3*initial_direction+i);imagesc(reshape(aligned(:,i),numThetan,nTau+1));
% end
