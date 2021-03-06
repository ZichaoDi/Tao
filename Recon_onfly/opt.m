%%% Simulate XRF of a given fixed object with rotating predifined detector and beam
%%% Travelling of fluorescence photon is approximated as the area covered by solid angle
global low up penalty Tik lambda
global W0 current_n Wold
global Beta TempBeta beta_d
global err0 fiter nit maxiter
global ConstSub MU_XTM I0 DisR
global itertest ErrIter icycle maxOut
global COR initialize mtol xbox ybox L in_after
global LogScale Tol 
COR=0;
% slice=66;
do_setup;
more off;
icycle=0;
Define_Detector_Beam_Gaussian; %% provide the beam source and Detectorlet
DefineObject_Gaussian; % Produce W, MU_XTM
%%%----------------------------Initialize dependent variables
DecomposedElement=1;
if(strcmp(sample,'Seed'))
    beta_d=1e0;
elseif(strcmp(sample,'Rod'))
    beta_d=1e1;
else
    beta_d=1e2;
end
if(~synthetic & ~ReconAttenu)
    truncChannel=1*(DecomposedElement==0);
    if(DecomposedElement)
        XRF=XRF_decom;%./max(XRF_decom(:));
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
end
fprintf(1,'====== With Self Absorption, Transmission Detector Resolution is %d\n',nTau);
%%%----------------------------------------------------------------------
W0=W(:);
%%%================= First-order derivative regularization
penalty=0;
cmap=[min(W0),max(W0)];
Mt=-log(data_abs./max(data_abs(:)));
windowSize=2;
xrt1=zeros(size(Mt));
for i=1:numThetan, 
    xrt1(i,:)=medfilt1(Mt(i,:),windowSize);
end
Mt=xrt1(:);
Mt=Mt*beta_d;
clear y xrt1;

xrfData=XRF(:);
fctn_f=@(W)Forward_onfly(W,xrfData,Mt,M);
rng('default');
x0 = 0*10^(0)*rand(m(1)*m(2)*NumElement,1);
err0=norm(W0-x0);
NF = [0*N; 0*N; 0*N];
tic;
x=x0;
Wold=reshape(x,prod(m),NumElement);%zeros(prod(m),NumElement);
err=[];
if(~ReconAttenu && Alternate)
    % initialize=1;
    %%%%%%%==============================================================
    if(initialize)
        mtol=prod(m);
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        L=sparse(zeros(numThetan*(nTau+1),mtol));%zeros(numThetan*(nTau+1)*prod(m),1);
        in_after=cell(numThetan*(nTau+1),mtol);
    end
    [ConstSub,fold]=feval(fctn_f,x);
    initialize=0;
    res0=fold;
    res=[];
end
fctn=@(W)sfun_linear_joint(W,xrfData,Mt,MU_e,L,NumElement,m,nTau);
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
        maxOut=1;
    end
    fprintf(1, 'cycle       alpha         residual      error      sub-residual\n');
    while(icycle<=maxOut)
        if(maxOut>1)
            figure(12);
            x_temp=reshape(x,m(1),m(2),NumElement);
            for ie=1:NumElement;
                subplot(maxOut+1,NumElement,(icycle-1)*NumElement+ie),
                imagesc(x_temp(:,:,ie))
                set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
            end
            clear x_temp
            drawnow;
        end
        maxiter=60;
        if(bounds)
            [x,f,g,ierror] = tnbc (x,fctn,low,up); % algo='TNbc';
        else
            [x,f,g,ierror] = tn (x,fctn);
        end
        [ConstSub,fold]=feval(fctn_f,x);
        % [x,fold,ConstSub, alpha]=lin3(x-Wold(:),Wold(:),fold,1,fctn_f,ConstSub);
        alpha=1;
        Wold=x;
        err(icycle)=norm(W0-x);
        res(icycle)=fold;
        fprintf(1,'%4i      %.2e       %.3e     %.3e    %.3e\n', icycle,alpha, fold, err(icycle),f);
        err_Joint=err;
        figure(10);h1=plot(0:icycle,[res0,res],'r.-');drawnow;hold on;
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
     save(['xstar_',sample,'_',num2str(N(1)),'_',num2str(numThetan),'_',num2str(TempBeta),'_',num2str(Beta),'_',num2str(beta_d),'_',num2str(numChannel),'-',num2str(nit),frame,num2str(slice),'_',num2str(icycle),'.mat'],'x_admm');
end
%%%%%%%%%%%%%%%%%%%%======================================================
if(~ReconAttenu)
    t=toc;
    fprintf('time elapsed = %f \n',t);
end
if(plotResult)
    doplot(1,xstar,W0,x0);
end
