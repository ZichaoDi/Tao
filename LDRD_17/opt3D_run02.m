global N_delta maxiter W0 sinoS
global xiter fiter ErrIter

ReconAttenu = 1; % 0: Recover W; 1: Recover miu
slice=1;
tau_rate=1;
ang_rate=1;
do_setup;
N_delta=numThetan;
rng('default');
res=3;
x_res=[];
aligned=[];
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,1);
initial_direction=size(d0,2);
deltaStar=zeros(N_delta,1);
W0=[deltaStar;W(:)];
maxiter=1000;
for element_ind=1:NumElement
recon0=zeros(N^2,nslice);
recon1=zeros(N^2+numThetan,nslice);
    for slice=1:nslice

        load(['Run02_',num2str(slice)]);
        XRF_decom=-log(data(1:tau_rate:end,1:ang_rate:end)'./max(max(data)));
        % setup_paunesku;
        % XRF_decom=squeeze(XRF_decom3D(:,slice,:,:));
        errW=zeros(size(d0,2),1);
        f_global=zeros(size(d0,2),1);
        NF = [0*N; 0*N; 0*N];
        Mt=XRF_decom(:,:,element_ind)';
        Mt=Mt./max(Mt(:));%%==normalize data;
        sinoS=Mt;
        Lmap=sparse(L);
        low=[-nTau/2*ones(N_delta,1);zeros(prod(m)*NumElement,1)];
        up=[nTau/2*ones(N_delta,1);inf*ones(prod(m)*NumElement,1)];
        W0=W(:);
        err0=norm(W0);
        fctn=@(x)sfun_radon(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
        [x,f,g,ierror] = tnbc (zeros(N^2*NumElement,1),fctn,low(N_delta+1:end),up(N_delta+1:end)); % algo='TNbc';
        W0=[deltaStar;W(:)];
        for res_step=1:initial_direction
            delta=(cos(theta)-1).*d0(1,res_step)+sin(theta).*d0(2,res_step);
            x0= [delta;0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
            err0=norm(W0-x0);
            fctn_COR=@(x)sfun_cor_dr(x,full(Mt'),sparse(Lmap));% on attenuation coefficients miu;
            bounds=1;
            if(bounds)
                [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
                errW(res_step)=norm(W0(N_delta+1:end)-xCOR(N_delta+1:end));
                f_global(res_step)=f;
                [~,~,alignedSignal]=feval(fctn_COR,xCOR);
                err_sino(res_step)=norm(alignedSignal'-sinoS);
                %% ============================================
            else
                [xCOR,f,g,ierror] = tn (x0,fctn_COR);
            end
            x_res=[x_res,xCOR];
            aligned=[aligned,alignedSignal(:)];
        end
        recon1(:,slice)=xCOR;recon0(:,slice)=x;
    % save(['result/olga/recon_reduced',num2str(slice),'_',num2str(N),'_',num2str(numThetan),'.mat'],'x','xCOR');
        
    end
    save(['result/run02/recon_reduced',num2str(element_ind),'_',num2str(ang_rate),'_',num2str(tau_rate),'.mat'],'recon1','recon0');
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
