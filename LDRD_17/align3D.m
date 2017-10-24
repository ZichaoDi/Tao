global N_delta maxiter W0 sinoS
global xiter fiter ErrIter
ReconAttenu = 1; % 0: Recover W; 1: Recover miu
slice=1;
% do_setup;
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
        % Mt=Mt-min(Mt(:));
        sinoS=Mt;
        Lmap=sparse(L);
        low=[-inf.*ones(N_delta,1);zeros(prod(m)*NumElement,1)];
        up=[inf.*ones(N_delta,1);inf*ones(prod(m)*NumElement,1)];
        for res_step=1:initial_direction
            delta=(cos(theta)-1).*d0(1,res_step)+sin(theta).*d0(2,res_step);
            x0= [delta;0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
            err0=norm(x0-x0);
            fctn_COR=@(x)sfun_cor_dr(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
            % load aligned86_5full;
            % shift=recon(1:N_delta,slice);
            x0(1:N_delta)=shift(:,slice);
            [~,~,aligned_xtm]=feval(fctn_COR,x0);
        end
        noAligned(:,:,slice,element_ind)=XRF_decom(:,:,element_ind);
        aligned(:,:,slice,element_ind)=aligned_xtm;
    end
end
save(['result/paunesku_align/aligned_mc.mat'],'aligned','noAligned');
