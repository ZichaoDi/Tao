
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
%%%----------------------------------------------------------------------

x0_XRF=x0;
x0_XTM=sum(reshape(x0,m(1),m(2),NumElement).*repmat(MUe,[m(1),m(2),1]),3);
x0_XTM=x0_XTM(:);
xinitial=x0_XRF;
err0_XTM=norm(x0_XTM(:)-MU_XTM(:));
err0_XRF=norm(x0_XRF-WS(:));
errTolOld=err0_XRF;
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
    NF = [0*N; 0*N; 0*N];
    Joint=(-1)^(OuterIter+1); % 1: XRF; -1: XTM; 0: Joint inversion
    if(Joint==-1)
        fprintf('============================ Transmission Reconstruction\n')
        maxiter=10;
        W0=MU_XTM(:);
        err0=norm(x_XTM(:)-MU_XTM(:));
        fctn=@(MU)sfun_XTM_com(DisR,MU,I0,L,thetan,m,nTau);
        [x_XTM,f,g,ierror] = tnbc (x_XTM(:),fctn,low_XTM,up_XTM);
    elseif(Joint==1)
        fprintf('============================ Fluorescence Reconstruction\n')
        maxiter=50;
        W0=W(:);
        err0=norm(x_XRF-WS(:));
        fctn=@(W)sfun_Tensor4_diff(W,x_XTM,XRF,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau);
        [x_XRF,f,g,ierror] = tnbc (x_XRF,fctn,low_XRF,up_XRF);
        x_XTM=sum(reshape(x_XRF,m(1),m(2),NumElement).*repmat(MUe,[m(1),m(2),1]),3);
    end
    
    %%%========================================================================
    if(Joint==1)
        Ntot(2)=Ntot(2)+NF(2)+NF(3);
    elseif(Joint==-1)
        Ntot(1)=Ntot(1)+NF(2)+NF(3);
    end
    
    OuterIter=OuterIter+1;
    if(Joint==1)
        errTol=norm(x_XRF-WS(:));
        if(abs(errTol-errTolOld)<eps)
            break;
        end
        errTolOld=errTol;
    elseif(Joint==-1)
        errTol=norm(x_XTM-MU_XTM(:));
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
