global N_delta maxiter W0 sinoS
global xiter fiter ErrIter
global alpha lambda
global reg_str

initialize=0;
if(initialize)
    sample='chip';
    ReconAttenu = 0; % 0: Recover W; 1: Recover miu
    ind_cr=1;
    synthetic=0;
    do_setup;
else
    load work
    setup_chip;
end
N_delta=numThetan;
rng('default');
res=3;
d0=[0 -res -res 0     res res res 0   -res;...
    0  0   -res -res -res 0   res res res];
d0=d0(:,1);
initial_direction=size(d0,2);
deltaStar=zeros(N_delta,1);
W0=[deltaStar;repmat(W(:),nslice,1)];
maxiter=100;
aligned3D=zeros(numThetan*(nTau+1),NumElement,nslice);
recon=sparse(N^2*NumElement+N_delta,nslice);
NF = [0*N; 0*N; 0*N];
rng('default')
Mt= permute(data_h,[3 2 1]);
    prj=0*Mt;
    for ele=1:nslice
        a=Mt(:,:,ele);alpha(ele)=std(a(:));
        for n=1:numThetan
            Mt(n,:,ele)=map1D(Mt(n,:,ele),[0,1]);
        end
        prj(:,:,ele)=Mt(:,:,ele);
        % Mt(:,:,ele)=Mt(:,:,ele)./max(max(Mt(:,:,ele)));
    end
    MtPad=GaussianPadded(Mt.^(1/1),0*max(nTau+1,nchannel),0);
    % save('pydata.mat','prj','thetan');
    % return;
    alpha=(alpha./sum(alpha));
    sinoS=reshape(Mt,numThetan*(nTau+1),nslice);
    Lmap=sparse(L);
%     for res_step=1:initial_direction
%         delta=repmat(d0(:,res_step),n_delta/2,1);
%         if(initialize)
%             x0 = [(cos(theta)-1)*delta(1)+sin(theta)*delta(2);0*10^(0)*rand(m(1)*m(2)*nslice,1)];
%         else
%             x0=xCOR;%[zeros(numThetan,1);x_base(:)];
%         end
%          err0=norm(W0-x0);
%          fctn_COR=@(x)sfun_shift_mc(x,MtPad,Lmap);% on attenuation coefficients miu;
%          low=[-inf/2*ones(N_delta,1);zeros(prod(m)*nslice,1)];
%          up=[inf/2*ones(N_delta,1);inf*ones(prod(m)*nslice,1)];
%          bounds=1;
%         reg_str={'TV'};
%         lambda=0e-5;
%         lam=-6;%[-15:2:-4];%, 0.2 0.4 0.8 1.6];
%         for i=1:length(lam)
%             lambda=0*10^lam(i);
%             [~,~,alignedSignal]=feval(fctn_COR,x0);
%              if(bounds)
%                  [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
% %                 errW(res_step)=norm(W0(n_delta+1:end)-xCOR(n_delta+1:end));
% %                 f_global(res_step)=f;
%                  [~,~,alignedSignal]=feval(fctn_COR,xCOR);
% %             else
% %                 [xCOR,f,g,ierror] = tn (x0,fctn_COR);
%              end
%              return;
%         end
%     end
% % prj=reshape(alignedSignal,numThetan,nTau+1,nslice);
% % if(initialize)
    reg_str={'TV'};
    lambda=1e-4;
    shift=zeros(numThetan,1);
    Mt0=reshape(Mt,numThetan*(nTau+1),nchannel);
    x_base=zeros(N^2,nchannel);
    maxOuter=5;
    maxiter=100;
    for iter=1:maxOuter
        demo_l2_TV;
        x_base=x_twist;
        % for ele=1:nchannel
        %     W0=zeros(N^2,1);
        %     err0=norm(W0-x_base(:,ele));
        %     fctn=@(x)sfun_radon(x,reshape(Mt0(:,ele),numThetan,nTau+1),Lmap);
        %     x_base(:,ele)=tnbc(x_base(:,ele),fctn,zeros(N^2,1),inf*ones(N^2,1));
        % end
        fctn_p=@(p)sfun_p(p,x_base,Mt,Lmap);
        shift=tn(shift,fctn_p);
        [~,~,Mt0]=feval(fctn_p,shift);
    end
% end
% save(['result/',sample,'/Recon_shift_mc_alter',num2str(numThetan),'_',num2str(nchannel),'_',num2str(nTau+1),'.mat'],'x_base','Mt0','shift');

nrow=3;
figure,
xc=reshape(x_base,[N,N,nchannel]);%reshape(xCOR(N_delta+1:end,:),[N,N,nslice,initial_direction]);
alignedSignal=Mt0;
for ele=1:nslice
    subplot(nrow,nslice,ele);
    imagesc(Mt(:,:,ele));
    set(gca,'xtick',[],'ytick',[]);
    if(ele==1);xlabel('\tau');ylabel('original: \theta');end
    subplot(nrow,nslice,nslice+ele);
    align=reshape(alignedSignal,numThetan*(nTau+1),nslice);
    imagesc(reshape(align(:,ele),numThetan,nTau+1));
    set(gca,'xtick',[],'ytick',[]);

    if(ele==1);xlabel('\tau');ylabel('aligned: \theta');end
    subplot(nrow,nslice,nslice*2+ele);
    imagesc(xc(:,:,ele));
    if(ele==1);ylabel('reconstruction');end
end

