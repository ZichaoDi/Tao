global TakeLog
global NF N ptest gv
global low up
close all;
TakeLog=0;
more off;
XRF_XTM_Gaussian;
% save('data/3x3_90XRF_FULL.mat','W','XRF','MU_e','M','NumElement','numChannel','Ltol','GlobalInd','thetan','m','nTau');
% load 3x3_90XRF_FULL;
% return;
%%%----------------------------------------------------------------------
W0=W;
%   fctn=@(W)sfun_XRF(W,XRF,MU,M,NumElement,numChannel,Ltol,GlobalInd,thetan,m,nTau);
if(TakeLog)
fctn=@(W)sfun_XRF_EM(W,XRF,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,thetan,m,nTau);
else
fctn=@(W)sfun_XRF_full3(W,XRF,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,thetan,m,nTau);
end
rng('default');
%%%================================================================
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
%        gs=foo(fctn,x1);%+rand(prod(m)*size(M,1),1)
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
        %         err(jj,i)=norm(fH(i)-f)/(1+norm(f));
        %         fprintf('%4.8e           %4.8e           %4.2e\n',f,fH(i),h(i));
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
 return;
end
% % % % %
  x0=W(:)+1*10^(0)*rand(prod(m)*size(M,1),1);
  xinitial=x0;
%   [f,g]=feval(fctn,x0);
%   t=cputime-e
% %   return;
%      foo(fctn,x0);
%      return;
% load x_new
% x0=x_new;
% % %    x0=W(:)+ones(prod(m)*size(M,1),1);
%  options = optimset('Algorithm','interior-point','DerivativeCheck','off','Diagnostics','off','GradObj','on','Display','iter');%,'AlwaysHonorConstraints','none','TolCon',1e-10,'TolX',1e-16,'TolFun',1e-15
% % [x,fval] = fminunc(fctn,x0,options);
% % Aeq=zeros(prod(m),length(x0)); 
% % for it=1:prod(m)
% %     for ij=1:NumElemen t
% %             Aeq(it,it+(ij-1)*prod(m))=1;
% %     end
% % end
% % beq=ones(prod(m),1); 
%  [xstar,fval] = fmincon(fctn,x0,[],[],[],[],zeros(size(x0)),[],[],options);
%  return;
% bounds=1;
N = 3;
NF = [0*N; 0*N; 0*N];
e=cputime;
figureObject(reshape(x0,m(1),m(2),NumElement),Z,m,NumElement,MU_e,1);
low=0*ones(size(x0));
up=1e6*ones(size(x0));
[xstar,f,g,ierror] = tnbc (x0,fctn,low,up);
figureObject(reshape(xstar,m(1),m(2),NumElement),Z,m,NumElement,MU_e,2);
%%%%%%%%%%%%%%%%%%%%=======================================================
figure(24);
ws=W(:);
for i=1:NumElement
    subplot(3,NumElement,i);
    plot(1:prod(m),x0(prod(m)*i-prod(m)+1:prod(m)*i),'ro',1:prod(m),xstar(prod(m)*i-prod(m)+1:prod(m)*i),'bs',1:prod(m),ws(prod(m)*i-prod(m)+1:prod(m)*i),'g*')
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
  err=norm(xstar-W0(:))/norm(xinitial-W0(:));
t=cputime-e;
fprintf('residule is %d, time elapsed is %d \n',err,t);
