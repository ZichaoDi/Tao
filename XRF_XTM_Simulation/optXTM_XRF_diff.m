
global maxiter NF N current_n Ntot
global ptest gv Joint VarInd W0 err0
close all;
maxiter=10;
LogScale=1;
XRF_XTM_Gaussian;
%%%----------------------------------------------------------------------
W0=W(:);
rng('default');
Wtest=W;
ws=Wtest(:);
x0_XRF=Wtest(:)+1*10^(0)*rand(prod(m)*size(M,1),1);%10*ones(sizEDe(ws));%
x0_XTM=MU_XTM(:)+1*10^(0)*rand(prod(m),1);
xinitial=x0_XRF;
err0=xinitial-ws;
err0_XTM=x0_XTM-MU_XTM(:);
N=m(1);
current_n=N(1);
e=cputime;
figureObject(reshape(x0_XRF,m(1),m(2),NumElement),Z,m,NumElement,MU_e,1);
low_XRF=0*ones(size(x0_XRF));
up_XRF=1e6*ones(size(x0_XRF));
low_XTM=0*ones(size(x0_XTM));
up_XTM=1e6*ones(size(x0_XTM));
%%%===== Initialization
OuterIter=0;
errTol=1;
x_XRF=x0_XRF;
x_XTM=x0_XTM;
Ntot=zeros(1,2);
VarInd=1:length(W0);%[10:18];
%%%===================================================================
MaxiOuter=2^2+1; %% choose odd number to make sure the outer iteration stops at XRF
while ( OuterIter<=MaxiOuter);
    NF = [0*N; 0*N; 0*N];
    Joint=(-1)^(OuterIter+1); % 1: XRF; -1: XTM; 0: Joint inversion
    if(Joint==-1)
        fctn=@(MU)sfun_XTM_com(DisR,MU,I0,Ltol,thetan,m,nTau);
        [x_XTM,f,g,ierror] = tnbc (x_XTM,fctn,low_XTM,up_XTM);
    elseif(Joint==1)
        fctn=@(W)sfun_XRF_full3(W,XRF,x_XTM,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,LocalInd,L_after,thetan,m,nTau);
        [x_XRF,f,g,ierror] = tnbc (x_XRF,fctn,low_XRF,up_XRF);
    end
    
    %%%========================================================================
    if(Joint==1)
        Ntot(2)=Ntot(2)+NF(2)+NF(3);
    elseif(Joint==-1)
        Ntot(1)=Ntot(1)+NF(2)+NF(3);
    end
    
    OuterIter=OuterIter+1;
    if(Joint==1)
        errTol=norm(x_XRF-ws)/norm(err0);
        if(errTol<1e-10)
            break;
        end
    elseif(Joint==-1)
        errTol=norm(x_XTM-MU_XTM(:))/norm(err0_XTM);
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
xstar=x_XRF;

for i=1:NumElement
    err(i)=norm(xstar(9*i-8:9*i)-ws(9*i-8:9*i))/norm(xinitial(9*i-8:9*i)-ws(9*i-8:9*i));
end
figure('name','Elemental Residule')
semilogy(1:NumElement,err,'r.-');
t=cputime-e;
fprintf('Time elapsed is %f, residule is %d\n',t,errTol);
%%%%%%%%%%%%%%%%%%%%=======================================================
figureObject(reshape(xstar,m(1),m(2),NumElement),Z,m,NumElement,MU_e,2);
%%%%%%%%%%%%%%%%%%%%=======================================================
figure(24);
for i=1:NumElement
    subplot(3,NumElement,i);
    plot(1:prod(m),xinitial(prod(m)*i-prod(m)+1:prod(m)*i),'ro',1:prod(m),xstar(prod(m)*i-prod(m)+1:prod(m)*i),'bs',1:prod(m),ws(prod(m)*i-prod(m)+1:prod(m)*i),'g*')
    xlim([0 prod(m)]);
    if(i==1)
        legend('initial','final','optimal','font',16)
        ylabel('solution','fontsize',12)
    end
    title(['Element ',num2str(i)],'fontsize',12);
end
for i=1:NumElement; subplot(3,NumElement,i+NumElement);plot(1:prod(m),ptest(prod(m)*i-prod(m)+1:prod(m)*i),'r.','MarkerSize',12);
    xlim([0 prod(m)]);
    if(i==1)
        ylabel('Projected Direction','fontsize',12)
    end
end
for i=1:NumElement; subplot(3,NumElement,i+2*NumElement);plot(1:prod(m),gv(prod(m)*i-prod(m)+1:prod(m)*i),'r.','MarkerSize',12);
    xlim([0 prod(m)]);
    if(i==1)
        ylabel('Projected Gradient','fontsize',12)
    end
end
%%%%%%%%%%%%%%%%%%%%=======================================================
