for slice= 86;%1:nslice
    NF = [0*N; 0*N; 0*N];
    Mt=XRF_raw(:,:,:,slice);
    for ele=1:NumElement
        a=Mt(:,:,ele);alpha(ele)=std(a(:));
        Mt(:,:,ele)=Mt(:,:,ele)./max(max(Mt(:,:,ele)));
    end
    alpha=alpha./sum(alpha);
    % Mt=Mt./max(Mt(:))+abs(min(Mt(:)));%%==normalize data;
    sinoS=reshape(Mt,numThetan*(nTau+1),NumElement);
    Lmap=sparse(L);
    for res_step=1:initial_direction
        delta=repmat(d0(:,res_step),n_delta/2,1);
        x0 = [(cos(theta)-1)*delta(1)+sin(theta)*delta(2);0*10^(0)*rand(m(1)*m(2)*NumElement,1)];
        err0=norm(W0-x0);
        fctn_COR=@(x)sfun_shift_mc(x,Mt,Lmap);% on attenuation coefficients miu;
        low=[-nTau/2*ones(N_delta,1);zeros(prod(m)*NumElement,1)];
        up=[nTau/2*ones(N_delta,1);inf*ones(prod(m)*NumElement,1)];
        bounds=1;
        % load aligned86_5_paunesku;
        % x0(1:numThetan)=shift;
        reg_str={'L2'};
        lam=-15;%[-15:2:-4];%, 0.2 0.4 0.8 1.6];
        recon=sparse(N^2*NumElement+N_delta,length(lam));
        for i=1:length(lam)
            lambda=0e-6;%10^lam(i);
            [~,~,alignedSignal]=feval(fctn_COR,x0);
            if(bounds)
                [xCOR,f,g,ierror] = tnbc (x0,fctn_COR,low,up); % algo='TNbc';
                errW(res_step)=norm(W0(n_delta+1:end)-xCOR(n_delta+1:end));
                f_global(res_step)=f;
                [~,~,alignedSignal]=feval(fctn_COR,xCOR);
                %% ============================================
                % x1=zeros(N^2,NumElement);
                % for ele=1:NumElement
                %     for t=1:numThetan
                %         alignedSignal(t:numThetan:end,ele)=medfilt1(alignedSignal(t:numThetan:end,ele),6);
                %     end

                %     fctn=@(x)sfun_radon(x,full(alignedSignal(:,ele)),sparse(Lmap));% on attenuation coefficients miu;
                %     W0=zeros(N^2,1);
                %     [x,f,g,ierror] = tnbc (W0,fctn,zeros(N^2,1),inf*ones(N^2,1)); % algo='TNbc';
                %     x1(:,ele)=x;
                % end
                %% ============================================
            else
                [xCOR,f,g,ierror] = tn (x0,fctn_COR);
            end
            recon(:,i)=xCOR;
            aligned3D(:,:,i)=alignedSignal;
        end
    end
end
% save(['result/',sample,'/Recon_shift_mc_weighted_',reg_str{1},'.mat'],'lam','recon','aligned3D');

nrow=3;
figure,
xc=reshape(xCOR(N_delta+1:end,:),[N,N,NumElement,initial_direction]);
for ele=1:NumElement
    subplot(nrow,NumElement,ele);
    imagesc(XRF_raw(:,:,ele,slice));
    title(Element{Z(ele)});
set(gca,'xtick',[],'ytick',[]);
    if(ele==1);xlabel('\tau');ylabel('original: \theta');end
    subplot(nrow,NumElement,NumElement+ele);
    align=reshape(alignedSignal,numThetan*(nTau+1),NumElement);
    imagesc(reshape(align(:,ele),numThetan,nTau+1));
set(gca,'xtick',[],'ytick',[]);

    if(ele==1);xlabel('\tau');ylabel('aligned: \theta');end
    subplot(nrow,NumElement,NumElement*2+ele);
    imagesc(xc(:,:,ele));
    if(ele==1);ylabel('reconstruction');end
end
