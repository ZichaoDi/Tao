
global NF N current_n ptest gv 
global low up Tik penalty gama  
global W0 LogScale maxiter err0 Joint

close all;
more off;
PlotObject=1;
plotElement=1;
LogScale=1;
maxiter=10000;
XRF_XTM_Gaussian;
%%%----------------------------------------------------------------------
W0=W(:);
Joint=1; % 1: XRF; -1: XTM; 0: Joint inversion
%%%============== Rescale MU_e to make unity contribution
DiscreteScale=0;
penalty=0;

if(DiscreteScale )
    Mag=-order(MU_e(:,1,1));
    gama=-log(MU_e(:,1,1))./MU_e(:,1,1);%5e2*ones(size(MU_e(:,1,1)));%10.^(Mag);%(max_MU-min_MU)./(MU_e(:,1,1)-min_MU);%1./MU_e(:,1,1);%
else
    gama=ones(size(MU_e(:,1,1)));
end

if(penalty)
    Tik=spalloc(length(W(:))-1,length(W(:)),2*(length(W(:))-1)); %% First-order derivative regularization
    for i=1:length(W(:))-1
        Tik(i,i)=-1; Tik(i,i+1)=1;
    end
    
    % Tik=spalloc(length(W(:))-2,length(W(:)),3*(length(W(:))-2)); %% First-order derivative regularization
    %     for i=1:size(Tik,1)
    %         Tik(i,i)=-1;Tik(i,i+1)=2; Tik(i,i+2)=-1;
    %     end
end

if(Joint==-1)
    fctn=@(W)sfun_XTM(W,DisR,MU_e,I0,Ltol,thetan,m,nTau,NumElement);
elseif(Joint==0)
    fctn=@(W)sfun_XRF_XTM(W,XRF,DisR,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,LocalInd,L_after,thetan,m,nTau,I0);
%      fctn=@(W)sfun_XRF_XTM_Jacobian(W,XRF,DisR,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,LocalInd,L_after,thetan,m,nTau,I0);
else
     fctn=@(W)sfun_XRF_full2(W,XRF,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,LocalInd,L_after,thetan,m,nTau);
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
% foo(fctn,x0);
% return;
N=m(1);
current_n=N(1);
NF = [0*N; 0*N; 0*N];
e=cputime;
figureObject(reshape(x0,m(1),m(2),NumElement),Z,m,NumElement,MU_e,1);
low=0*ones(size(x0));
up=1e6*ones(size(x0));
%%%========================================================================
[xstar,f,g,ierror] = tnbc (x0,fctn,low,up);
% xstar_Joint1=xstar;
% save xstar_Joint1 xstar_Joint1;
%%%========================================================================
% load xstar_Joint1
% xstar=xstar_Joint1;
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
% errTol=(norm(xstar-ws)^2/length(xstar))^(1/2)/(max(ws)-min(ws)); %% NRMSE 
if(DiscreteScale)
    AbsErr=norm(xtemp(:)-W(:));
    IniErr=norm(err0_1);
    errOri=norm(xtemp(:)-W(:))/norm(err0_1);
    fprintf('Time elapsed is %f, residule is %d, original residule is %d\n',t,errTol,errOri);
else
    fprintf('Time elapsed is %f, residule is %d\n',t,errTol);
end
%%%%%%%%%%%%%%%%%%%%=======================================================
% load xstar_Joint1
% load xstar_diff1
% for iplot=2:3
%     if(iplot==2)
% xstar=xstar_Joint1;
%     elseif(iplot==3)
%         xstar=xstar_diff1;
%     end
% figureObject(reshape(xstar,m(1),m(2),NumElement),Z,m,NumElement,MU_e,iplot);
% end
%%%%%%%%%%%%%%%%%%%%=======================================================
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
%     plot(1:prod(m),ptest(prod(m)*i-prod(m)+1:prod(m)*i),'r.','MarkerSize',12);
%     xlim([0 prod(m)]);
%     if(i==1)
%         ylabel('Projected Direction','fontsize',12)
%     end
end
for i=1:NumElement; subplot(3,NumElement,i+2*NumElement);plot(1:prod(m),sort(gv(prod(m)*i-prod(m)+1:prod(m)*i)),'r.','MarkerSize',12);
    xlim([0 prod(m)]);
    if(i==1)
        ylabel('Projected Gradient','fontsize',12)
    end
end
% %%%%%%%%%%%%%%%%%%%%=======================================================
