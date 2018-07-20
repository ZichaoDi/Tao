global N_delta maxiter W0 sinoS
global xiter fiter ErrIter
ReconAttenu = 1; % 0: Recover W; 1: Recover miu
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
for element_ind=1;% :length(slice_tot)
    recon=sparse(N^2,nslice);
    for slice=1:nslice
        % setup_miller;
        NF = [0*N; 0*N; 0*N];
        Mt=squeeze(XRF_decom3D(:,slice,:))';%XRF_decom(:,:,element_ind)';
        sinoS=Mt;
        Lmap=sparse(L);
        low=[-inf.*ones(N_delta,1);zeros(prod(m)*NumElement,1)];
        up=[inf.*ones(N_delta,1);inf*ones(prod(m)*NumElement,1)];
        for res_step=1:initial_direction
            delta=(cos(theta)-1).*d0(1,res_step)+sin(theta).*d0(2,res_step);
            x0= [delta;0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
            err0=norm(x0-x0);
            fctn_COR=@(x)sfun_cor_dr(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
            % load aligned30_2_miller

            x0(1:N_delta)=shift_mean;
            [~,~,aligned_xtm]=feval(fctn_COR,x0);
             % for i=1:numThetan
             %     aligned_xtm(i,:)=map1D(aligned_xtm(i,:),[0,1]);
             % end
             % aligned_xtm=aligned_xtm./max(aligned_xtm(:));
             aligned(:,:,slice,element_ind)=aligned_xtm;
            fctn=@(x)sfun_radon(x,full(aligned_xtm),sparse(Lmap));% on attenuation coefficients miu;

            bounds=1;
            [x,f,g,ierror] = tnbc (x0(N_delta+1:end),fctn,low(N_delta+1:end),up(N_delta+1:end)); % algo='TNbc';
        end
        recon(:,slice)=x;
    end
    % save(['result/miller/recon_aligned_raw',num2str(element_ind),'.mat'],'recon','aligned');
    % save(['result/miller/recon_aligned_norm',num2str(element_ind),'.mat'],'recon','aligned');
end
