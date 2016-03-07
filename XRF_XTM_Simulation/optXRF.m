global low up penalty
global W0 current_n Wold
global SigMa_XTM SigMa_XRF Beta TempBeta
global err0 fiter nit maxiter
global ConstSub InTens OutTens MU_XTM
global itertest ErrIter
more off;
%%%----------------------------Initialize dependent variables
% do_setup;
level=1;
current_n = N(level);
W= W_level{level};
DecomposedElement=0;
if(synthetic)
    XRF=double(xrf_level{level});
else
    if(DecomposedElement)
        XRF=double(xrf_level_decom{level});
        DetChannel=DetChannel_decom;
        numChannel=numChannel_decom;
        M=M_decom*4e1;
    else
        XRF=double(xrf_level_raw{level});
        DetChannel=DetChannel_raw;
        numChannel=numChannel_raw;
        M=M_raw;
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
if(NoSelfAbsorption)
    fprintf(1,'====== No Self Absorption, Transmission Detector Resolution is %d\n',nTau);
else
    fprintf(1,'====== With Self Absorption, Transmission Detector Resolution is %d\n',nTau);
end
%%%----------------------------------------------------------------------
W0=W(:);
cmap=[min(W0),max(W0)];
if(TakeLog)
    fctn=@(W)sfun_XRF_EM(W,XRF,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,thetan,m,nTau);
elseif(Alternate)
    fctn_f=@(W)Calculate_Attenuation(W,NumElement,L,GlobalInd,SelfInd,m,nTau,XRF,M);
    if(Joint==0)
        TempBeta=1; Beta=0;
        if(linear_S)
            Mrep=repmat(reshape(M,[1,1,1 NumElement numChannel]),[numThetan,nTau+1,mtol,1,1]);
            fctn=@(DW)sfun_linearized(DW,XRF,Mrep,NumElement,L,m,nTau);
        else
            Mrep_P=repmat(reshape(M,[1,1,1 NumElement numChannel]),[numThetan,nTau+1,mtol,1,1]);
            fctn_half=@(W)sfun_half_linear(W,XRF,M,NumElement,L,GlobalInd,SelfInd,m,nTau);
            L_rep=reshape(L,numThetan,nTau+1,mtol);
            fctn=@(W)sfun_full_linear(W,XRF,NumElement,m,nTau);
        end

    elseif(Joint==1)
        TempBeta=1; Beta=1e7;
        Mrep_P=repmat(reshape(M,[1,1 NumElement numChannel]),[nTau+1,mtol,1,1]);
        fctn=@(W)sfun_Joint_admm(W,XRF,DisR,MU_e,Mrep_P,NumElement,L,GlobalInd,SelfInd,m,nTau);
        % fctn_f=@(W)func_Tensor_Joint(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0);
    end
else
    if(Joint==0)
        TempBeta=1; Beta=0;
        fctn=@(W)sfun_Tensor(W,XRF,M,NumElement,L,GlobalInd,SelfInd,m,nTau);
    elseif(Joint==1)
        TempBeta=1; Beta=1e9;
        fctn=@(W)sfun_Tensor_Joint(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,m,nTau,I0);
    elseif(Joint==-1)
        TempBeta=0; Beta=1;
        if(ReconAttenu)
            fctn=@(MU)sfun_XTM_tensor(DisR,MU,I0,L,m,nTau);% on attenuation coefficients miu;
        else
            fctn=@(W)sfun_XTM(W,DisR,MU_e,I0,L,m,nTau,NumElement);%on W
        end
    end
end
rng('default');
x0 = 0*W(:)+1*10^(1)*rand(prod(m)*NumElement,1);
xinitial=x0;
err0=norm(W0-xinitial);
NF = [0*N; 0*N; 0*N];
tic;
maxOut=3;
icycle=1;
x=x0;
Wold=reshape(x,prod(m),NumElement);%zeros(prod(m),NumElement);
err=[];
errDW=[];
errAtt=[];
Wtemp=[];
outCycle=1;
if(~ReconAttenu)
    [InTens, OutTens, AttenuM, DW,fold]=Calculate_Attenuation(Wold,NumElement,L,GlobalInd,SelfInd,m,nTau,XRF,M);
    temp_v=L_rep.*InTens;
    ConstSub=bsxfun(@times,bsxfun(@times,temp_v,OutTens),Mrep_P); % part 3
    clear DW temp_v;
    % [InStar, OutStar, AttenuStar, DWstar,fstar]=Calculate_Attenuation(W,NumElement,L,GlobalInd,SelfInd,m,nTau,XRF,M);
    % errAtt0=norm(AttenuM(:)-AttenuStar(:));
end
%%%==================================================
low=0*ones(size(x0));
up=1e6*ones(size(x0));
x_admm=[];
x_admm(:,1)=x0;
err_obj=[];
if(Alternate & linear_S==0)
    maxOut=1;
    x=x0;
    % RealVec=[];
    % for ib=1:size(RealBeam,2),RealVec(ib)=RealBeam(2,ib)+(nTau+1)*(RealBeam(1,ib)-1);end
    fctn0=@(W)sfun_Tensor(W,XRF,M,NumElement,L,GlobalInd,SelfInd,m,nTau);
    fprintf(1, 'cycle       alpha         residual      error      sub-residual\n');
    while(icycle<=maxOut)
    %%%%% =================== Attenuation Matrix at beam energy
        matrixVersion=0;
        if(matrixVersion)
            MU_XTM=sum(reshape(Wold,m(1),m(2),NumElement).*repmat(MUe,[m(1),m(2),1]),3).*XTMscale;
            d=XRF(:);
            if(strcmp(grad_type,'full-linear'))
                Cmatrix=reshape(L,numThetan,nTau+1,mtol);
                B=bsxfun(@times,bsxfun(@times,Cmatrix.*InTens,OutTens),reshape(M,1,1,1,NumElement,numChannel));
                B=permute(B,[1,2,5,3,4]);
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
            % RealVec=find(d~=0);
            % B=B(RealVec,:);
            % d=d(RealVec);
            % x=solver_ART(B,x,d,2*size(B,1));%B(RealVec,:)\d(RealVec); algo='ART';
            % [x,flag,relres,iter]=pcg(B'*B,B'*d,1e-6,size(B,2)); algo='pcg';
            x=B\d; algo='backslash';
            fctn_linear=@(x)linear(x,B'*B,B'*d);%B(RealVec,:),d(RealVec));
            f=norm(B*x-d);
        end
        maxiter=40;
        if(bounds)
            [x,f,g,ierror] = tnbc (x,fctn,low,up); algo='TNbc';
            % options = optimoptions('fmincon','Display','iter','GradObj','on','MaxIter',maxiter,'Algorithm','interior-point');% ,'Algorithm','Trust-Region-Reflective'
            % x=fmincon(fctn,x,[],[],[],[],low,up,[],options); algo='fmincon';
        else
            [x,f,g,ierror] = tn (x,fctn);
        end
        [x,fold,InTens, OutTens,AttenuM, alpha]=lin3(x-Wold(:),Wold(:),fold,1,fctn_f,InTens, OutTens, AttenuM);
        temp_v=L_rep.*InTens;
        ConstSub=bsxfun(@times,bsxfun(@times,temp_v,OutTens),Mrep_P); % part 3
        clear DW temp_v;
        x_admm(:,icycle+1)=x;
        Wold=x;
        % AttenuOld=AttenuStar;
        err(icycle)=norm(W0-x);
        % errAtt(icycle)=norm(AttenuM(:)-AttenuStar(:));
        fprintf(1,'%4i      %.2e       %.3e     %.3e    %.3e\n', icycle,alpha, fold, err(icycle),f);
        if(Joint==0)
            err_XRF=err;
            plot_str='ko-';
            figure(10);h1=semilogy(0:icycle,[err0,err_XRF],plot_str);
            text(icycle,err_XRF(icycle),['f = ',num2str(fold,'%2.2e')]);
            drawnow;hold on;%
            % hAtt=semilogy(0:icycle,[errAtt0,errAtt],'b*-');
        elseif(Joint==1)
            err_Joint=err;
            figure(10);h1=semilogy(0:icycle,[err0,err_Joint],'r.-');drawnow;hold on;
            text(icycle,err_Joint(icycle),['f = ',num2str(fold,'%2.2e')]);
            drawnow;hold on;%
            % hAtt=semilogy(0:icycle,[errAtt0,errAtt],'b*-');
        end
        if(icycle>1| icycle==maxOut | fold<=1e-5)
            if(icycle==maxOut | norm(x-W0)<1e-6 );%| err(icycle)>=err(icycle-1))
                xstar=x;
                break;
            end
        end
        icycle=icycle+1;
    end
    % h_plot=[h_plot,h1];
    save(['xstar_',sample,'_',num2str(N(1)),'_',num2str(numThetan),'_',num2str(numel(num2str(TempBeta))),'_',num2str(numel(num2str(Beta))),'_',num2str(numChannel),'alternate.mat'],'xstar','x0','W0');
    % legend([h1;hAtt],'error-W','error-Attenuation');

elseif(Alternate==0)
    maxiter=100;
    if(bounds)
        [xstar,f,g,ierror]=tnbc(x,fctn,low,up);
    else
        [xstar,f,g,ierror]=tn(x,fctn);
    end
    % save(['xstar_',sample,'_',num2str(N(1)),'_',num2str(numThetan),'_',num2str(numel(num2str(TempBeta))),'_',num2str(numel(num2str(Beta))),'_full.mat'],'xstar');
    if(Joint==0)
        ErrIter_XRF=ErrIter;
        figure(10); h3=semilogy(linspace(0,maxOut,length(ErrIter_XRF)), ErrIter_XRF,'ms-');hold on;
        % legend([h;h1;h3],'error-W','Alternate-tn','Full-tn')
    elseif(Joint==1)
        ErrIter_Joint=ErrIter;
        figure(10); h4=semilogy(linspace(0,maxOut,length(ErrIter_Joint)), ErrIter_Joint,'g*-');hold on;
    end
end
%%%%%%%%%%%%%%%%%%%%=======================================================
if(~ReconAttenu)
    t=toc;
    fprintf('time elapsed = %f \n',t);
end
%%%%%%%%%%%%%%%%%%%%=======================================================
if(plotResult)
    doplot(1,xstar,W_level,xinitial);
end
%%%======================= Derivative Test
TestDirection=0;
if(TestDirection)
    h=linspace(-0.5,12,50);
    hKnot=[1:3:12];%[1:5];%[-0.01:0.01:0.05];% 
    ntot=length(hKnot);
    Np=1;
    for jj=1:Np
        plus=5*10^(-1)*rand(prod(m)*size(M,1),1);%ones(prod(m)*size(M,1),1);%
        fs=[];
        gss=[];
        for t=1:length(hKnot)
        x1=W0(:)+hKnot(t)*plus;
        [fs(t),gs]=feval(fctn,x1);
        gss(t)=gs'*plus;
        end
        fH=[];
        fs1=[];
        for i=1:length(h)
            x0=W0(:)+h(i)*plus;
            [fH(i),gH] =feval(fctn,x0);
            for t=1:length(hKnot)
            fs1(i,t)=fs(t)+gss(t)*(h(i)-hKnot(t));
            end
        end
        figure,
        subplot(Np,ntot+1,(ntot+1)*jj-ntot);plot(plus,'r.-');
        if(jj==1)
            xlabel('h');ylabel('w_{r}','Interpreter','Tex');
        end
        for t=1:ntot
            subplot(Np,ntot+1,(ntot+1)*jj-(ntot-t)),plot(h,fH,'r-',h,fs1(:,t),'b-');title(['h=',num2str(hKnot(t))]);
            if(jj==1 & t==1)
                xlabel('h');legend('f(w*+hw_{r})','f(w_1)+(h-h_1)w_{r}^{T} \nabla f(w_1)','Interpreter','latex');
            end
        end
    end
end
%%%=============================== Linear system starts
if(linear_S)
    xstar=squeeze(DW(RealBeam(1,1),RealBeam(2,1),:));
    [low,up]=set_bounds(1,1);
    C=reshape(L,numThetan*(nTau+1),prod(m));
    BT=kron(M',C);
    xsT=solver_ART(BT,0*xstar,XRF(:),size(RealBeam,2));%BT\XRF(:);
    while (icycle<=maxOut)
        for im=1:size(RealBeam,2)
            RB=RealBeam(:,im);
            Mrep=repmat(reshape(M,[1,1,1 NumElement numChannel]),[numThetan,nTau+1,mtol,1,1]);
            nT=RB(2);
            nAngle=RB(1);
            fprintf(1,'Angle = %f,  Beamlet = %d\n',thetan(nAngle),nT);
            Dstar=DWstar(nAngle,nT,:);
            d=squeeze(XRF(nAngle,nT,:));
            Cmatrix=reshape(L(nAngle,nT,:,:),1,prod(m));
            B=kron(M',Cmatrix);
            xtemp=solver_ART(B,xstar,d,1);% B\d;
            % fctn=@(DW)sfun_linearized_single(DW,squeeze(XRF(nAngle,nT,:)),squeeze(Mrep(nAngle,nT,:,:,:)),NumElement,squeeze(L(nAngle,nT,:,:)),m);
            % [xtemp,f,g,ierror] = tnbc (xstar,fctn,low,up);
            % fprintf(1,'==================================================\n');
            InterSectCount=0*ones(prod(m),1);
            Index=GlobalInd{nAngle,nT};
            InterSectCount(sub2ind(m,Index(:,2),Index(:,1)))=1;
            xstar = xtemp.*(xtemp>=0);% xstar+InterSectCount.*xtemp;
            W_update=squeeze(reshape(xstar,1,1,mtol,NumElement)./AttenuM(nAngle,nT,:,:));
            Wold=(1-repmat(InterSectCount,1,NumElement)).*Wold+repmat(InterSectCount,1,NumElement).*W_update;%repmat(1./IntersectCount,1,NumElement).*% 
            [InTens, OutTens, AttenuM, DW,f]=Calculate_Attenuation(Wold,NumElement,L,GlobalInd,SelfInd,m,nTau,XRF,M);
            Wtemp{icycle,nAngle,nT}=Wold;
            errDW(im)=norm(xstar-squeeze(Dstar));
            err(im)=norm(Wold(:)-W(:));
            errAtt(im)=norm(AttenuM(:)-AttenuStar(:));
        end
        W1=repmat(reshape(xsT,1,1,mtol,NumElement),[numThetan,nTau+1,1,1])./AttenuStar;
        figure(9);
        if(icycle==1)
            h=semilogy(0+im*(icycle-1):im+im*(icycle-1),[err0,err],'ro-',1+im*(icycle-1):im+im*(icycle-1),errDW,'b*-',1+im*(icycle-1):im+im*(icycle-1),errAtt,'g.-');
        else
            h=semilogy(1+im*(icycle-1):im+im*(icycle-1),err,'ro-',1+im*(icycle-1):im+im*(icycle-1),errDW,'b*-',1+im*(icycle-1):im+im*(icycle-1),errAtt,'g.-');
        end
        legend('err-W','err-DW','err-Att');
        hold on;
        icycle=icycle+1;
    end
        figure(10);h=semilogy(linspace(0,icycle,length(err)+1),[err0,err],'ro-');
        hold on;
    xstar=Wold(:);
end
