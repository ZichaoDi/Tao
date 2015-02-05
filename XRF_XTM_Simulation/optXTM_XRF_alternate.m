
global maxiter NF N low up current_n
global ptest gv
global penalty gama
global W0 err0
global SigMa_XTM SigMa_XRF

%%%----------------------------Initialize dependent variables
more on;
plotResult=1;
do_setup;
level=find(N==current_n);
W= W_level{level};
XRF=xrf_level{level};
DisR=xtm_level{level};
L=L_level{level};
GlobalInd=GI_level{level};
SelfInd=SI_level{level};
m=m_level(level,:);
nTau=nTau_level(level);
SigMa_XTM=SigmaT{level};
SigMa_XRF=SigmaR{level};

%%%----------------------------------------------------------------------
W0=W(:);
%%%============== Rescale MU_e to make unity contribution
DiscreteScale=0;
penalty=0;
xinitial=x0;
e=cputime;
figureObject(reshape(x0,m(1),m(2),NumElement),Z,m,NumElement,MU_e,1);
low=0*ones(size(x0));
up=1e6*ones(size(x0));
%%%===== Initialization
OuterIter=0;
errTol=1;
x=x0;
Ntot=zeros(1,2);
%%%===================================================================
while (OuterIter <= 15);
    NF = [0*N; 0*N; 0*N];
    Joint=(-1)^OuterIter; % 0: XRF; -1: XTM; 1: Joint inversion
    err0=norm(x-W0(:));
    if(Joint==-1)
        maxiter=10;
        fprintf('============================ Transmission Reconstruction\n');
        fctn=@(W)sfun_XTM(W,DisR,MU_e,I0,L,thetan,m,nTau,NumElement);
    elseif(Joint==1)
        maxiter=50;
        fprintf('============================ Fluorescence Reconstruction\n');
        fctn=@(W)sfun_Tensor4(W,XRF,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau);
    end
    errTolOld=norm(x-W0);
    %%%========================================================================
    [x,f,g,ierror] = tnbc (x,fctn,low,up);
    OuterIter=OuterIter+1;
    errTol=norm(x-W0);
    
    if(Joint==1)
        Ntot(2)=Ntot(2)+NF(2)+NF(3);
        if(abs(errTol-errTolOld)<1e-6)
            break;
        end
    elseif(Joint==-1)
        Ntot(1)=Ntot(1)+NF(2)+NF(3);
    end
    figure(333)
    hold on; drawnow;
    if(Joint==1)
        marker='s';
        color='r';
    elseif(Joint==-1)
        marker='*';
        color='b';
    end
    plot(OuterIter,errTol,'Marker',marker,'color',color,'LineStyle','-')
end
fprintf('Total number of function evaluations, XTM = %d, XRF = %d \n',Ntot(1),Ntot(2));
%%%====================================================== Report Result
xstar=x;
if(Joint==-1 & DiscreteScale)
    xtemp=reshape(xstar,m(1),m(2),NumElement);
    err0_1=reshape(err0,m(1),m(2),NumElement);
    for i=1:NumElement
        xtemp(:,:,i)=xtemp(:,:,i)*gama(i);
        err0_1(:,:,i)=err0_1(:,:,i)*gama(i);
    end
    err0_1=err0_1(:);
end
for i=1:NumElement
    err(i)=norm(xstar(9*i-8:9*i)-W0(9*i-8:9*i))/norm(xinitial(9*i-8:9*i)-W0(9*i-8:9*i));
end
figure('name','Elemental Residule')
semilogy(1:NumElement,err,'r.-');
t=cputime-e;
if(DiscreteScale)
    errOri=norm(xtemp(:)-W(:))/norm(err0_1);
    fprintf('Time elapsed is %f, residule is %d, original residule is %d\n',t,errTol,errOri);
else
    fprintf('Time elapsed is %f, residule is %d\n',t,errTol);
end
%%%%%%%%%%%%%%%%%%%%=======================================================
figureObject(reshape(xstar,m(1),m(2),NumElement),Z,m,NumElement,MU_e,2);
%%%%%%%%%%%%%%%%%%%%=======================================================
if(plotResult)
    figure(24);
    for i=1:NumElement
        subplot(3,NumElement,i);
        
        errCom=reshape(xstar(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-W0(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(errCom);colormap gray
        if(i==1)
            ylabel('Final Soluction','fontsize',12)
        end
        title(['Element ',num2str(i)],'fontsize',12);
    end
    
    for i=1:NumElement
        subplot(3,NumElement,i+NumElement);
        
        errCom=reshape(W0(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));
        imagesc(errCom);colormap gray
        if(i==1)
            ylabel('True Soluction','fontsize',12)
        end
    end
    for i=1:NumElement; subplot(3,NumElement,i+2*NumElement);
        plot(1:prod(m),sort(xinitial(prod(m)*i-prod(m)+1:prod(m)*i)),'ro',1:prod(m),sort(xstar(prod(m)*i-prod(m)+1:prod(m)*i)),'bs',1:prod(m),sort(W0(prod(m)*i-prod(m)+1:prod(m)*i)),'g*')
        xlim([0 prod(m)]);
        if(i==1)
            legend('initial','final','optimal','font',16)
            ylabel('solution','fontsize',12)
        end
    end
    % for i=1:NumElement; subplot(3,NumElement,i+2*NumElement);plot(1:prod(m),sort(gv(prod(m)*i-prod(m)+1:prod(m)*i)),'r.','MarkerSize',12);
    %     xlim([0 prod(m)]);
    %     if(i==1)
    %         ylabel('Projected Gradient','fontsize',12)
    %     end
    % end
end
%%%%%%%%%%%%%%%%%%%%=======================================================
