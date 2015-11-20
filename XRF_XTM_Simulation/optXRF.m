global low up penalty
global W0 current_n Wold
global SigMa_XTM SigMa_XRF Beta TempBeta
global fctn_f err0 fiter nit maxiter
%%%----------------------------Initialize dependent variables
% do_setup;
level=1;
current_n = N(level);
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
if(NoSelfAbsorption)
    fprintf(1,'====== No Self Absorption, Transmission Detector Resolution is %d\n',nTau);
else
    fprintf(1,'====== With Self Absorption, Transmission Detector Resolution is %d\n',nTau);
end
%%%----------------------------------------------------------------------
W0=W(:);
if(TakeLog)
    fctn=@(W)sfun_XRF_EM(W,XRF,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,thetan,m,nTau);
elseif(Alternate)
    fctn=@(W)sfun_Tensor_Patrick(W,XRF,M,NumElement,L,GlobalInd,SelfInd,m,nTau);
    fctn1=@(W)sfun_Tensor(W,XRF,M,NumElement,L,GlobalInd,SelfInd,m,nTau);
else
    fctn=@(W)sfun_Tensor(W,XRF,M,NumElement,L,GlobalInd,SelfInd,m,nTau);
end
rng('default');
x0=W(:)+1*10^(0)*rand(prod(m)*size(M,1),1);
xinitial=x0;
NF = [0*N; 0*N; 0*N];
e=cputime;
low=0*ones(size(x0));
up=1e6*ones(size(x0));
maxiter=20;
maxOut=1;
icycle=1;
x=x0;
Wold=x0;
while(icycle<=maxOut)
    [x,f,g,ierror] = tnbc (x,fctn,low,up);
    Wold=x;
    if(icycle==maxOut | norm(x-W0)<1e-6)
        xstar=x;
        break;
    end
    icycle=icycle+1;
end
%%%%%%%%%%%%%%%%%%%%=======================================================
err=norm(xstar-W0(:))/norm(xinitial-W0(:));
t=cputime-e;
fprintf('residule is %d, time elapsed is %d \n',err,t);
save(['xs',num2str(N(1)),'_',num2str(numThetan),'Golosio',num2str(TempBeta),'_',num2str(Beta),'_BI.mat'],'xstar','x0');
%%%%%%%%%%%%%%%%%%%%=======================================================
if(plotResult)
    figure('name','Elemental Residule');
     clims=[0 max(W0)];
    for i=1:NumElement
        subplot(4,NumElement,i);
        
        errCom=reshape(xinitial(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-W0(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(errCom,clims);colormap jet
        if(i==1)
            ylabel('Initial Guess','fontsize',12)
        end
        title(Element{Z(i)},'fontsize',12);
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
    end
    
    for i=1:NumElement
        subplot(4,NumElement,i+1*NumElement);
        
        errCom=reshape(xstar(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-W0(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(errCom,clims);colormap jet
        if(i==1)
            ylabel('Final Soluction','fontsize',12)
        end
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
    end
    
    for i=1:NumElement
        subplot(4,NumElement,i+2*NumElement);
        
        errCom=reshape(W0(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));
        imagesc(errCom,clims);colormap jet
        if(i==1)
            ylabel('True Soluction','fontsize',12)
        end
        if(i==NumElement)
            hp4 = get(subplot(4,NumElement,i+2*NumElement),'Position');
            colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(3)]);
        end
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
    end
    for i=1:NumElement; subplot(4,NumElement,i+3*NumElement);
        plot(1:prod(m),sort(xinitial(prod(m)*i-prod(m)+1:prod(m)*i)),'ro',1:prod(m),sort(xstar(prod(m)*i-prod(m)+1:prod(m)*i)),'bs',1:prod(m),sort(W0(prod(m)*i-prod(m)+1:prod(m)*i)),'g*')
        xlim([0 prod(m)]);
        if(i==1)
            hleg=legend('initial','final','optimal','FontSize',6, 'Box', 'off');
            set(hleg,'units','pixels');
            lp=get(hleg,'outerposition');
            set(hleg,'Location','NorthWest', 'Box', 'off','outerposition',[lp(1),10,lp(3),20]);
            ylabel('solution','fontsize',12)
        end
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
    end
    
    % for i=1:NumElement; subplot(3,NumElement,i+2*NumElement);plot(1:prod(m),sort(gv(prod(m)*i-prod(m)+1:prod(m)*i)),'r.','MarkerSize',12);
    %     xlim([0 prod(m)]);
    %     if(i==1)
    %         ylabel('Projected Gradient','fontsize',12)
    %     end
    % end
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
