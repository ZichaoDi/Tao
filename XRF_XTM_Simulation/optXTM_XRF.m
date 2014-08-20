
global NF N ptest gv
global low up 
global min_MU max_MU DiscreteScale
close all;
more off;
% SimulateXTM;
XRF_XTM_Gaussian;
DiscreteScale=0;
%%%----------------------------------------------------------------------
W0=W;
Joint=1; % 1: XRF; -1: XTM; 0: Joint inversion
if(Joint==-1)
    fctn=@(W)sfun_XTM(W,DisR,MU_e,I0,Ltol,thetan,m,nTau,NumElement);
elseif(Joint==0)
    fctn=@(W)sfun_XRF_XTM(W,XRF,DisR,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,thetan,m,nTau,I0);
else
    fctn=@(W)sfun_XRF_full3(W,XRF,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,LocalInd,L_after,thetan,m,nTau);
end
rng('default');
Wtest=W;
for i=1:NumElement
    Wtest(:,:,i)=Wtest(:,:,i)/(-log(MU_e(i,1,1)))*MU_e(i,1,1);%*(MU_e(i,1,1)-min_MU)/(max_MU-min_MU);
end
   x0=Wtest(:)*1e-2+1*10^(1)*rand(prod(m)*size(M,1),1);
%    x0=ones(size(x0));
  xinitial=x0;
  
%    options = optimset('Algorithm','interior-point','DerivativeCheck','off','Diagnostics','off','GradObj','on','Display','iter');%,'AlwaysHonorConstraints','none','TolCon',1e-10,'TolX',1e-16,'TolFun',1e-15
%  [xstar,fval] = fmincon(fctn,x0,[],[],[],[],zeros(size(x0)),ones(size(x0)),[],options);
% return;
%      foo(fctn,x0);
%      return;
% load x_new
% x0=x_new;
N=m(1);
NF = [0*N; 0*N; 0*N];
e=cputime;
figureObject(reshape(x0,m(1),m(2),NumElement),Z,m,NumElement,MU_e,1);
low=0*ones(size(x0));
up=1e6*ones(size(x0));
[xstar,f,g,ierror] = tnbc (x0,fctn,low,up);
figure(99);plot(xstar)
if(Joint==-1 & DiscreteScale)
xtemp=reshape(xstar,m(1),m(2),NumElement);
for i=1:NumElement
    xtemp(:,:,i)=xtemp(:,:,i)*(-log(MU_e(i,1,1)))/MU_e(i,1,1);%*(MU_e(i,1,1)-min_MU)/(max_MU-min_MU);
end
xstar=xtemp(:);
end
%%%====================================================== Report Result
Wstar=W0(:);
for i=1:NumElement
err(i)=norm(xstar(9*i-8:9*i)-Wstar(9*i-8:9*i))/norm(xinitial(9*i-8:9*i)-Wstar(9*i-8:9*i));
end
figure('name','Elemental Residule')
semilogy(1:NumElement,err,'r.-');
t=cputime-e;
errTol=norm(xstar-Wstar)/norm(xinitial-Wstar);
fprintf('Time elapsed is %f, residule is %d\n',t,errTol);
%%%%%%%%%%%%%%%%%%%%=======================================================
figureObject(reshape(xstar,m(1),m(2),NumElement),Z,m,NumElement,MU_e,2);
%%%%%%%%%%%%%%%%%%%%=======================================================
figure(24);
ws=W(:);
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
