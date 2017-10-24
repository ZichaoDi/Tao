global N_delta maxiter W0 sinoS
global xiter fiter ErrIter

ReconAttenu = 1; % 0: Recover W; 1: Recover miu
slice=1;
do_setup;
rng('default');
W0=W(:);
maxiter=150;
for element_ind=1:length(slice_tot)
    recon=sparse(N^2,nslice);
    for slice=1:nslice
        setup_miller;
        NF = [0*N; 0*N; 0*N];
        Mt=XRF_decom(:,:,element_ind)';
        Mt=Mt./max(Mt(:));
        Mt=Mt-min(Mt(:));
        sinoS=Mt;
        Lmap=sparse(L);
        low=zeros(prod(m)*NumElement,1);
        up=inf*ones(prod(m)*NumElement,1);
        x0= [0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
        err0=norm(x0-x0);
        fctn=@(x)sfun_radon(x,full(sinoS'),sparse(Lmap));% on attenuation coefficients miu;
        bounds=1;
        [x,f,g,ierror] = tnbc (x0,fctn,low,up); % algo='TNbc';
        recon(:,slice)=x;
    end
    save(['result/miller/recon0',num2str(element_ind),'.mat'],'recon');
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
