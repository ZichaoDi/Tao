global low up penalty Tik lambda
global W0 current_n Wold
global Beta TempBeta beta_d
global err0 fiter nit maxiter
global ConstSub InTens OutTens MU_XTM I0 s_a DisR
global itertest ErrIter icycle maxOut
global pert
more off;
%%%----------------------------Initialize dependent variables
if(strcmp(sample,'Seed'))
    beta_d=1e0;
elseif(strcmp(sample,'Rod'))
    beta_d=1e1;
end

if(n_level==1)
    if(~synthetic & ~ReconAttenu)
        truncChannel=1*(DecomposedElement==0);
        if(DecomposedElement)
            XRF=XRF_decom;  
            DetChannel=DetChannel_decom;
            numChannel=numChannel_decom;
            M=M_decom;
        else
            XRF=XRF_raw;
            DetChannel=DetChannel_raw;
            numChannel=numChannel_raw;
            M=M_raw;
            if(truncChannel)
                XRF=XRF(:,truncInd);
                DetChannel=DetChannel_raw(truncInd);
                numChannel=length(truncInd);
                M=M_raw(:,truncInd);
            end
        end
        % s_a=0;
        % if(s_a)
        %     DisR=squeeze(data_sa)';
        % else
        %     DisR=squeeze(data_ds)';  
        % end
    end
        % DisR=squeeze((data_h(3,:,:)-data_h(2,:,:))./(data_h(1,:,:)-data_h(2,:,:)))';
else
    current_n = N(level);
    W= W_level{level};
    if(synthetic)
        XRF=double(xrf_level{level});
    else
        if(DecomposedElement)
            XRF=double(xrf_level_decom{level});  
            DetChannel=DetChannel_decom;
            numChannel=numChannel_decom;
            M=M_decom;
        else
            XRF=double(xrf_level_raw{level});
            DetChannel=DetChannel_raw;
            numChannel=numChannel_raw;
            M=M_raw;
            if(truncChannel)
                XRF=XRF(:,:,truncInd);
                DetChannel=DetChannel_raw(truncInd);
                numChannel=length(truncInd);
                M=M_raw(:,truncInd);
            end
        end
    end
    drift = 0;
    if(drift)
        DisR=xtm_level_pert{level};
    else
        s_a=0;
        if(s_a)
            DisR=squeeze(data_sa)';
        else
            DisR=squeeze(data_ds)';%xtm_level_true{level};%DisR_Simulated;%  
        end
    end
    L=L_level{level};
    GlobalInd=GI_level{level};
    SelfInd=SI_level{level};
    m=m_level(level,:);
    nTau=nTau_level(level);
end
fprintf(1,'====== With Self Absorption, Transmission Detector Resolution is %d\n',nTau);
%%%----------------------------------------------------------------------
W0=W(:);
%%%================= First-order derivative regularization
penalty=0;
if(penalty)
    lambda=1e-4;
    Tik=delsq(numgrid('S',N(1)+2)); 
end
cmap=[min(W0),max(W0)];
if(Alternate)
    if(Joint==0)
        TempBeta=1; Beta=0;
        fctn_half=@(W)sfun_half_linear(W,XRF,M,NumElement,L,GlobalInd,SelfInd,m,nTau);
        fctn=@(W)sfun_full_linear(W,XRF,NumElement,m,nTau);
    elseif(Joint==1)
        TempBeta=0; Beta=1;%=0.5;
        if(s_a)
            Mt=reshape(DisR',numThetan*(nTau+1),1);
        else
            Mt=-log(DisR'./I0)*beta_d; 
            Mt=Mt(:)-min(Mt(:));
        end
        xrfData=XRF(:);%/I0;%/reshape(repmat(I0,[1 1 numChannel]),numThetan*(nTau+1)*numChannel,1);
        fctn_f=@(W)Calculate_Attenuation_S1(W,NumElement,L,GlobalInd,SelfInd,m,nTau,xrfData,Mt,M);
        fctn=@(W)sfun_linear_joint(W,xrfData,Mt,MU_e,L,NumElement,m,nTau);
    end
else
    if(Joint==0)
        TempBeta=1; Beta=0;
        fctn=@(W)sfun_Tensor_Simplified(W,XRF,M,NumElement,L,GlobalInd,SelfInd,m,nTau);
    elseif(Joint==1)
        TempBeta=1; Beta=1e9;
        fctn=@(W)sfun_Tensor_Joint(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,m,nTau,I0);

    elseif(Joint==-1)
        TempBeta=1; Beta=0;
        if(s_a)
            Mt=reshape(DisR',numThetan*(nTau+1),1);
        else
            Mt=-log(DisR'./I0); 
            Mt=Mt(:)-min(Mt(:));
        end
        if(ReconAttenu)
            fctn=@(MU)sfun_XTM_MU(Mt,MU,I0,L,m,nTau);% on attenuation coefficients miu;
        else
            fctn=@(W)sfun_XTM(W,DisR,MU_e,I0,L,m,nTau,NumElement);%on W
        end
    end
end
rng('default');
x0 = 0*10^(0)*rand(m(1),m(2),NumElement);
x0=x0(:);
err0=norm(W0-x0);
NF = [0*N; 0*N; 0*N];
tic;
maxOut=3;
x=x0;
Wold=reshape(x,prod(m),NumElement);%zeros(prod(m),NumElement);
err=[];
errDW=[];
errAtt=[];
Wtemp=[];
outCycle=2;
if(~ReconAttenu && Alternate)
    [ConstSub, fold]=feval(fctn_f,Wold);
    res0=fold;
    res=[];
end
icycle = 1;
%%%==================================================
low=0*ones(size(x0));
up=inf*ones(size(x0));
x_admm=[];
x_admm(:,1)=x0;
err_obj=[];
if(Alternate && linear_S==0)
    if(TempBeta==0)
        maxOut=1;
    else
        maxOut=2;
    end
    x=x0;
    fprintf(1, 'cycle       alpha         residual      error      sub-residual\n');
    while(icycle<=maxOut)
    %%%%% =================== Attenuation Matrix at beam energy
        matrixVersion=0;
        if(matrixVersion)
            fctn0=@(W)sfun_Tensor(W,XRF,M,NumElement,L,GlobalInd,SelfInd,m,nTau);
            d=XRF(:);
            if(strcmp(grad_type,'full-linear'))
                B=permute(ConstSub,[1,2,5,3,4]);
                B=reshape(B,numThetan*(nTau+1)*numChannel,mtol*NumElement);
                plot_str='ko-';
            else
                [ftemp,gtemp,Jacob1,Jacob2,Jacob3]=feval(fctn0,Wold(:));
                dF=reshape(permute(reshape(Jacob1,numThetan*(nTau+1),mtol*NumElement,numChannel),[1 3 2]),numThetan*(nTau+1)*numChannel,mtol*NumElement);
                dA=reshape(permute(reshape(Jacob2,numThetan*(nTau+1),mtol*NumElement,numChannel),[1 3 2]),numThetan*(nTau+1)*numChannel,mtol*NumElement);
                dW=reshape(permute(reshape(Jacob3,numThetan*(nTau+1),mtol*NumElement,numChannel),[1 3 2]),numThetan*(nTau+1)*numChannel,mtol*NumElement);
                if(strcmp(grad_type,'half-linear'))
                    B=dW+dA;
                    plot_str='g*-';
                elseif(strcmp(grad_type,'original'))
                    B=dW+dA+dF;
                    plot_str='r.-';
                end
                d=d-dW*Wold(:)+B*Wold(:);
            end
            [x,x_debias,f,times,debias_start,mses,taus]= ...
                    SpaRSA(d,B,0);
            % x=B\d; algo='backslash';
            fctn_linear=@(x)linear(x,B'*B,B'*d);%B(RealVec,:),d(RealVec));
            f=norm(B*x-d);
        else
            if(maxOut>1)
                figure(12);
                x_temp=reshape(x,m(1),m(2),NumElement);
                for ie=1:NumElement; 
                    subplot(maxOut+1,NumElement,icycle*NumElement+ie),
                    imagesc(x_temp(:,:,ie))
                    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
                end
                clear x_temp
                drawnow;
            end
        maxiter=50;
        if(bounds)
            [x,f,g,ierror] = tnbc (x,fctn,low,up); % algo='TNbc';
            % options = optimoptions('fmincon','Display','iter','GradObj','on','MaxIter',maxiter,'Algorithm','interior-point');% ,'Algorithm','Trust-Region-Reflective'
            % [x, f] = fmincon(fctn,x,[],[],[],[],low,up,[],options); algo='fmincon';
        else
            [x,f,g,ierror] = tn (x,fctn);
        end
    end
        [x,fold,ConstSub, alpha]=lin3(x-Wold(:),Wold(:),fold,1,fctn_f,ConstSub);
        Wold=x;
        err(icycle)=norm(W0-x);
        res(icycle)=fold;
        fprintf(1,'%4i      %.2e       %.3e     %.3e    %.3e\n', icycle,alpha, fold, err(icycle),f);
        if(Joint==0)
            err_XRF=err;
            plot_str='ko-';
            figure(10);h1=semilogy(0:icycle,abs([res0,res]),plot_str);
            text(icycle,err_XRF(icycle),['f = ',num2str(err_XRF(icycle),'%2.2e')]);
            drawnow;
        elseif(Joint==1)
            err_Joint=err;
            figure(10);h1=plot(0:icycle,[res0,res],'r.-');drawnow;hold on;
        end
        x_admm(:,icycle+1)=x;
        if(icycle>1)
            if(icycle==maxOut || norm(x_admm(:,icycle)-x_admm(:,icycle-1))<1e-6 );%| err(icycle)>=err(icycle-1))
                xstar=x;
                break;
            end
        elseif(maxOut==1)
            xstar=x;
            break;
        end
        icycle=icycle+1;
    end
    save(['xstar_',sample,'_',num2str(N(1)),'_',num2str(numThetan),'_',num2str(TempBeta),'_',num2str(Beta),'_',num2str(numChannel),'-',num2str(nit),frame,'_linear_',num2str(lambda),'.mat'],'x_admm','W0');
elseif(Alternate==0)
    % pert=L*x0+log(DisR(:)./I0);
    pert = 0* DisR(:);
    for icycle =1;
    maxiter=100;
    if(bounds)
         [xstar,f,g,ierror]=tnbc(x,fctn,low,up);
         % options = optimoptions('fmincon','Display','iter','GradObj','on','MaxIter',maxiter,'Algorithm','interior-point');% ,'Algorithm','Trust-Region-Reflective'
         % [xstar, f] = fmincon(fctn,x,[],[],[],[],low,up,[],options); algo='fmincon';
    else
        [xstar,f,g,ierror]=tn(x,fctn);
    end
    % DR=-log(DisR'/I0);
    % pert=L*xstar+log(DisR(:)./I0);
end
    save(['xstar_',sample,'_',num2str(N(1)),'_',num2str(numThetan),'_',num2str(numel(num2str(TempBeta))),'_',num2str(numel(num2str(Beta))),'_',num2str(numChannel),'_full.mat'],'xstar','x0','W0');
    if(Joint==0)
        ErrIter_XRF=ErrIter;
        figure(10); h3=semilogy(linspace(0,maxOut,length(ErrIter_XRF)), ErrIter_XRF,'ms-');hold on;
    elseif(Joint==1)
        ErrIter_Joint=ErrIter;
        figure(10); h4=semilogy(linspace(0,maxOut,length(ErrIter_Joint)), ErrIter_Joint,'g*-');hold on;
    end
end
%%%%%%%%%%%%%%%%%%%%======================================================
if(~ReconAttenu)
    t=toc;
    fprintf('time elapsed = %f \n',t);
end
if(plotResult)
    doplot(1,xstar,W0,x0);
end
