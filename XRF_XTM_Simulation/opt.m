global low up penalty
global W0 current_n Wold
global SigMa_XTM SigMa_XRF Beta TempBeta
global err0 fiter nit maxiter
global ConstSub InTens OutTens MU_XTM
global itertest ErrIter icycle scaleMU maxOut
more off;
%%%----------------------------Initialize dependent variables
level=1;
current_n = N(level);
W= W_level{level};
DecomposedElement=1;
truncChannel=1*(DecomposedElement==0);
if(synthetic)
    XRF=double(xrf_level{level});
else
    if(DecomposedElement)
        XRF=double(xrf_level_decom{level});% XRF_Simulated_decom;% 
        % XRF=-log(XRF+1.001);
        DetChannel=DetChannel_decom;
        numChannel=numChannel_decom;
        M=M_decom;
        % M=-log(M_decom*1e0+1.001);
    else
        XRF=double(xrf_level_raw{level});
        DetChannel=DetChannel_raw;
        numChannel=numChannel_raw;
        M=M_raw*1e0;
        if(truncChannel)
            XRF=XRF(:,:,truncInd);
            DetChannel=DetChannel_raw(truncInd);
            numChannel=length(truncInd);
            M=M_raw(:,truncInd)*1e0;
        end
    end
end
DisR=xtm_level{level};
L=L_level{level};
GlobalInd=GI_level{level};
SelfInd=SI_level{level};
m=m_level(level,:);
nTau=nTau_level(level);
SigMa_XTM=SigmaT{level};
SigMa_XRF=SigmaR{level};
fprintf(1,'====== With Self Absorption, Transmission Detector Resolution is %d\n',nTau);
%%%----------------------------------------------------------------------
W0=W(:);
cmap=[min(W0),max(W0)];
if(Alternate)
    fctn_f=@(W)Calculate_Attenuation_Simplified(W,NumElement,L,GlobalInd,SelfInd,m,nTau,XRF,M);
    if(Joint==0)
        TempBeta=1; Beta=0;
        Mrep_P=repmat(reshape(M,[1,1,1 NumElement numChannel]),[numThetan,nTau+1,mtol,1,1]);
        fctn_half=@(W)sfun_half_linear(W,XRF,M,NumElement,L,GlobalInd,SelfInd,m,nTau);
        fctn=@(W)sfun_full_linear(W,XRF,NumElement,m,nTau);
    elseif(Joint==1)
        TempBeta=1; Beta=1e7;
        Mrep_P=repmat(reshape(M,[1,1 NumElement numChannel]),[nTau+1,mtol,1,1]);
        fctn=@(W)sfun_Joint_admm(W,XRF,DisR,MU_e,Mrep_P,NumElement,L,GlobalInd,SelfInd,m,nTau);
    end
else
    if(Joint==0)
        TempBeta=1; Beta=0;
        fctn=@(W)sfun_Tensor_Simplified(W,XRF,M,NumElement,L,GlobalInd,SelfInd,m,nTau);
    elseif(Joint==1)
        TempBeta=1; Beta=1e9;
        fctn=@(W)sfun_Tensor_Joint(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,m,nTau,I0);
    elseif(Joint==-1)
        TempBeta=0; Beta=1;
        if(ReconAttenu)
            fctn=@(MU)sfun_XTM_MU(DisR,MU,I0,L,m,nTau);% on attenuation coefficients miu;
        else
            fctn=@(W)sfun_XTM(W,DisR,MU_e,I0,L,m,nTau,NumElement);%on W
        end
    end
end
rng('default');
x0 =  0*W(:)+1*10^(0)*rand(prod(m)*NumElement,1);
xinitial=x0;
err0=norm(W0-xinitial);
NF = [0*N; 0*N; 0*N];
tic;
maxOut=3;
icycle=0;
x=x0;
Wold=reshape(x,prod(m),NumElement);%zeros(prod(m),NumElement);
err=[];
errDW=[];
errAtt=[];
Wtemp=[];
outCycle=1;
if(~ReconAttenu & Alternate)
    scaleMU=1e0;
    [ConstSub, ~, ~, ~, ~,fold]=Calculate_Attenuation_Simplified(0*Wold,NumElement,L,GlobalInd,SelfInd,m,nTau,XRF,M);
    res0=fold;
    res=[];
end
icycle = icycle + 1;
%%%==================================================
low=0*ones(size(x0));
up=1e6*ones(size(x0));
x_admm=[];
x_admm(:,1)=x0;
err_obj=[];
if(Alternate & linear_S==0)
    maxOut=1;
    x=x0;
    fctn0=@(W)sfun_Tensor(W,XRF,M,NumElement,L,GlobalInd,SelfInd,m,nTau);
    fprintf(1, 'cycle       alpha         residual      error      sub-residual\n');
    while(icycle<=maxOut)
    %%%%% =================== Attenuation Matrix at beam energy
        matrixVersion=0;
        if(matrixVersion)
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
        maxiter=40;
        if(bounds)
            clear Mrep_P
            % [x,f,g,ierror] = tnbc (x,fctn,low,up); algo='TNbc';
            options = optimoptions('fmincon','Display','iter','GradObj','on','MaxIter',maxiter,'Algorithm','interior-point');% ,'Algorithm','Trust-Region-Reflective'
            [x, f] = fmincon(fctn,x,[],[],[],[],low,up,[],options); algo='fmincon';
        else
            [x,f,g,ierror] = tn (x,fctn);
        end
    end
        x_pre=x;
        if(icycle==1)
            scaleMU=1e0;
        else
            scaleMU=1; %10^(6);
        end
        [x,fold,ConstSub, alpha]=lin3(x_pre-Wold(:),Wold(:),fold,1,fctn_f,ConstSub);
        x_admm(:,icycle+1)=x;
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
            figure(10);h1=semilogy(0:icycle,[err0,err_Joint],'r.-');drawnow;hold on;
            text(icycle,err_Joint(icycle),['f = ',num2str(fold,'%2.2e')]);
            drawnow;
        end
        if(icycle>1| icycle==maxOut | fold<=1e-5)
            if(icycle==maxOut | norm(x-W0)<1e-6 );%| err(icycle)>=err(icycle-1))
                xstar=x;
                break;
            end
        end
        icycle=icycle+1;
    end
    save(['xstar_',sample,'_',num2str(N(1)),'_',num2str(numThetan),'_',num2str(numel(num2str(TempBeta))),'_',num2str(numel(num2str(Beta))),'_',num2str(numChannel),'alternate.mat'],'xstar','x0','W0');
elseif(Alternate==0)
    maxiter=100;
    if(bounds)
        [xstar,f,g,ierror]=tnbc(x,fctn,low,up);
    else
        [xstar,f,g,ierror]=tn(x,fctn);
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
    doplot(1,xstar,W_level,xinitial);
end
