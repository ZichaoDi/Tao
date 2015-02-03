
global maxiter NF N current_n Ntot
global  Joint W0 err0
global SigMa_XTM SigMa_XRF
close all;
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
maxiter=10;
%%%----------------------------------------------------------------------
W0=W(:);
x0_XRF=W0+1*10^(0)*rand(prod(m)*size(M,1),1);%10*ones(sizEDe(W0));%
x0_XTM=ones(size(MU_XTM(:)+1*10^(0)*rand(prod(m),1)));
xinitial=x0_XRF;
err0=norm(xinitial-W0);
err0_XTM=x0_XTM-MU_XTM(:);
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
%%%===================================================================
MaxiOuter=2^6+1; %% choose odd number to make sure the outer iteration stops at XRF
while ( OuterIter<=MaxiOuter);
    Joint=(-1)^(OuterIter+1); % 1: XRF; -1: XTM; 0: Joint inversion
    if(Joint==-1)
        fctn=@(MU)sfun_XTM_com(DisR,MU,I0,L,thetan,m,nTau);
        [x_XTM,f,g,ierror] = tnbc (x_XTM,fctn,low_XTM,up_XTM);
    elseif(Joint==1)
        fctn=@(W)sfun_Tensor4(W,XRF,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau);
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
        errTol=norm(x_XRF-W0)/err0;
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
xstar_diff1=xstar;
% save xstar_diff1 xstar_diff1;
for i=1:NumElement
    err(i)=norm(xstar(9*i-8:9*i)-W0(9*i-8:9*i))/norm(xinitial(9*i-8:9*i)-W0(9*i-8:9*i));
end
figure('name','Elemental Residule')
semilogy(1:NumElement,err,'r.-');
t=cputime-e;
fprintf('Time elapsed is %f, residule is %d\n',t,errTol);
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
