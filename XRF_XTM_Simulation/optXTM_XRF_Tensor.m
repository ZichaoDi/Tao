
global NF N current_n ptest gv 
global low up penalty gama  
global W0 LogScale maxiter err0 Joint

close all;
more off;
PlotObject=1;
plotElement=1;
LogScale=1;
maxiter=10000;
XRF_XTM_Tensor;
%%%----------------------------------------------------------------------
W0=W(:);
Joint=0; % 1: XRF; -1: XTM; 0: Joint inversion
%%%============== Rescale MU_e to make unity contribution
DiscreteScale=0;
penalty=0;

if(DiscreteScale )
    Mag=-order(MU_e(:,1,1));
    gama=-log(MU_e(:,1,1))./MU_e(:,1,1);%5e2*ones(size(MU_e(:,1,1)));%10.^(Mag);%(max_MU-min_MU)./(MU_e(:,1,1)-min_MU);%1./MU_e(:,1,1);%
else
    gama=ones(size(MU_e(:,1,1)));
end

if(Joint==-1)
    fctn=@(W)sfun_XTM(W,DisR,MU_e,I0,Ltol,thetan,m,nTau,NumElement);
elseif(Joint==0)
    fctn=@(W)sfun_TensorJ1(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0);
else
     fctn=@(W)sfun_Tensor2(W,XRF,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau);
end
rng('default');
Wtest=W;
if(DiscreteScale)
    for i=1:NumElement
        Wtest(:,:,i)=Wtest(:,:,i)/gama(i);
        up=1e6*ones(size(W));
        up(:,:,i)=up(:,:,i)/gama(i);
    end
end
ws=Wtest(:);
x0=Wtest(:)+1*10^(-1)*rand(prod(m)*size(M,1),1);%10*ones(size(ws));%
xinitial=x0;
err0=xinitial-ws;
% tic;
% [f,g]=feval(fctn,W(:));
% toc;
% return;
foo(fctn,x0);
return;
N=m(1);
current_n=N(1);
NF = [0*N; 0*N; 0*N];
e=cputime;
figureObject(reshape(x0,m(1),m(2),NumElement),Z,m,NumElement,MU_e,1);
low=0*ones(size(x0));
up=1e6*ones(size(x0));
%%%========================================================================
%  options = optimset('Algorithm','interior-point','DerivativeCheck','off','Diagnostics','off','GradObj','on','Display','iter');%,'AlwaysHonorConstraints','none','TolCon',1e-10,'TolX',1e-16,'TolFun',1e-15
%  [xstar,fval] = fmincon(fctn,x0,[],[],[],[],zeros(size(x0)),[],[],options);
%  return;
 [xstar,f,g,ierror] = tnbc (x0,fctn,low,up);
if(Joint==-1 & DiscreteScale)
    xtemp=reshape(xstar,m(1),m(2),NumElement);
    err0_1=reshape(err0,m(1),m(2),NumElement);
    for i=1:NumElement
        xtemp(:,:,i)=xtemp(:,:,i)*gama(i);
        err0_1(:,:,i)=err0_1(:,:,i)*gama(i);
    end
    err0_1=err0_1(:);
end
%%%====================================================== Report Result

for i=1:NumElement
    err(i)=norm(xstar(9*i-8:9*i)-ws(9*i-8:9*i))/norm(xinitial(9*i-8:9*i)-ws(9*i-8:9*i));
end
% figure('name','Elemental Residule')
% semilogy(1:NumElement,err,'r.-');
t=cputime-e;
errTol=norm(xstar-ws)/norm(err0);
if(DiscreteScale)
    AbsErr=norm(xtemp(:)-W(:))
    IniErr=norm(err0_1)
    errOri=norm(xtemp(:)-W(:))/norm(err0_1);
    fprintf('Time elapsed is %f, residule is %d, original residule is %d\n',t,errTol,errOri);
else
    fprintf('Time elapsed is %f, residule is %d\n',t,errTol);
end

figure(24);
for i=1:NumElement
    subplot(3,NumElement,i);

errCom=reshape(xstar(prod(m)*i-prod(m)+1:prod(m)*i)-ws(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));
imagesc(errCom);colormap gray
    title(['Element ',num2str(i)],'fontsize',12);
end
for i=1:NumElement; subplot(3,NumElement,i+NumElement);
        plot(1:prod(m),sort(xinitial(prod(m)*i-prod(m)+1:prod(m)*i)),'ro',1:prod(m),sort(xstar(prod(m)*i-prod(m)+1:prod(m)*i)),'bs',1:prod(m),sort(ws(prod(m)*i-prod(m)+1:prod(m)*i)),'g*')
    xlim([0 prod(m)]);
    if(i==1)
        legend('initial','final','optimal','font',16)
        ylabel('solution','fontsize',12)
    end
end
for i=1:NumElement; subplot(3,NumElement,i+2*NumElement);plot(1:prod(m),sort(gv(prod(m)*i-prod(m)+1:prod(m)*i)),'r.','MarkerSize',12);
    xlim([0 prod(m)]);
    if(i==1)
        ylabel('Projected Gradient','fontsize',12)
    end
end
% %%%%%%%%%%%%%%%%%%%%=======================================================
