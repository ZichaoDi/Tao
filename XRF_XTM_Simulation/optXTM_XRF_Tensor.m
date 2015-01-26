
global NF N current_n ptest gv
global low up penalty gama
global W0 LogScale maxiter err0 Joint

close all;
more on;
PlotObject=1;
plotSpec = 0; % Do you want to see the spectra? If so plotSpec = 1
plotTravel=0; % If plot the intersection of beam with object
plotUnit=0;
plotElement=0;
LogScale=1;
maxiter=10000;
XRF_XTM_Tensor;
%%%----------------------------------------------------------------------
W0=W(:);
Joint=1; % 0: XRF; -1: XTM; 1: Joint inversion
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
    fctn=@(W)sfun_XTM(W,DisR,MU_e,I0,L,thetan,m,nTau,NumElement);
elseif(Joint==1)
%     fctn=@(W)sfun_Tensor_Joint_Jacobian(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0);
    fctn=@(W)sfun_Tensor_Joint(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0);
else
     fctn=@(W)sfun_Tensor4(W,XRF,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau);
    fctn1=@(W)sfun_AdiMat(W,XRF,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau);
end
Wtest=W;
if(DiscreteScale)
    for i=1:NumElement
        Wtest(:,:,i)=Wtest(:,:,i)/gama(i);
        up=1e6*ones(size(W));
        up(:,:,i)=up(:,:,i)/gama(i);
    end
end
ws=Wtest(:);
xinitial=x0;
err0=xinitial-ws;

%%%===================================================== Derivative Test
% gh=foo(fctn,x0);
% % % foo(fctn1,x0);
% % [f,g]=feval(fctn,x0);
% return;
% TolP=3;
% [f1,g1]=feval(fctn1,x0);
% gh=foo(fctn,x0);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% For loop version (previous)
% % disp('========= For Loop Version')
% % gFor=g;
% % XRF1=XRF;
% % optXTM_XRF;
% % %%%%%%%%%%%%==============================================
% f_AD_Tensor=abs(f-f1)/norm(f)
% % f_AD_For=norm(f1-fFor3)/norm(f)
% % f_Tensor_For=norm(f-fFor3)/norm(f)
% % return;
% errTensor_FD=norm(g-gh)/norm(g1)
% %  errFor_FD=norm(gFor3-gForh)/norm(g1)
% % errFD1_FD2=norm(g1-gFor3)/norm(g1)
% % return;
% errAD_Tensor=norm(g-g1)/norm(g1)
% % errFor_AD=norm(g1-gFor3)/norm(g1)
% % errFor_Tensor=norm(g-gFor3)/norm(g1)
% ee=abs(g-g1);
% % eeFor=abs(gFor-g1);
% figure('name','Derivative Difference'),
% subplot(TolP,1,1), hold on
% hl1= plot(1:prod(m)*NumElement,abs(g1-g)   ./ max(abs([g g1]),[],2),'r.-');
% hl2= plot(1:prod(m)*NumElement,abs(gh-g)./ max(abs([gh g]),[],2),'b.-');
% % hl3= plot(1:prod(m)*NumElement,abs(g-gFor3) ./ max(abs([gFor3 g]),[],2),'k.-');
% hl4= plot(1:prod(m)*NumElement,abs(g1-gh)  ./ max(abs([gh g1]),[],2),'g.-');
% % hl5= plot(1:prod(m)*NumElement,abs(gFor3-gForh)./max(abs([gForh gFor3]),[],2),'c.-');
% %  hl6= plot(1:prod(m)*NumElement,abs(g1-gForh)./max(abs([gForh g1]),[],2),'m.--');
% % hl7= plot(1:prod(m)*NumElement,abs(g1-gFor3)./max(abs([gForh3 g1]),[],2),'go--');
% 
% set(gca,'yscale','log'); set(gca,'ytick',10.^[-10:10])
% 
% axis tight
% % legend([hl1, hl2,hl3,hl4,hl5,hl6,hl7],'relerror(Tensor,AdiMat)','relerror(For2,AdiMat)',...
% %     'relerror(For,Tensor)','relerror(FD,AdiMat)','relerror(For,FDFor)','relerror(FDFor,AdiMat)','relerror(For3,AdiMat)')
% legend([hl1,hl2,hl4],'relerror(Tensor,AdiMat)',...
%     'relerror(FD,Tensor)','relerror(FD,AdiMat)')
% subplot(TolP,1,2),plot(1:prod(m)*NumElement,g,'r.-',1:prod(m)*NumElement,g1,'go-');%,1:prod(m)*NumElement,gFor,'b*-')
% legend('Tensor','AdiMat','For')
% hold on; for i=1:NumElement, line([prod(m)*i,prod(m)*i],[0,max(ee)],'LineStyle',':'); text(prod(m)*i,max(ee),num2str(i));end
% subplot(TolP,1,TolP),plot(sum(MU_e,3),'r.-');
% return;
% %%=======================================================================
N=m(1);
current_n=N(1);
NF = [0*N; 0*N; 0*N];
e=cputime;
low=0*ones(size(x0));
up=1e6*ones(size(x0));
%%%========================================================================
%  options = optimset('Algorithm','interior-point','Display','iter');%,'DerivativeCheck','off','Diagnostics','off','GradObj','off','Display','iter','AlwaysHonorConstraints','none','TolCon',1e-10,'TolX',1e-16,'TolFun',1e-15);%
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

% for i=1:NumElement
%     err(i)=norm(xstar(prod(m)*i-(prod(m)-1):prod(m)*i)-ws(prod(m)*i-(prod(m)-1):prod(m)*i))/norm(xinitial(prod(m)*i-(prod(m)-1):prod(m)*i)-ws(prod(m)*i-(prod(m)-1):prod(m)*i));
% end
% figure('name','Elemental Residule')
% semilogy(1:NumElement,err,'r.-');
t=cputime-e;
errTol=norm(xstar-ws)/norm(err0);
if(DiscreteScale)
    AbsErr=norm(xtemp(:)-W(:))
    IniErr=norm(err0_1)
    errOri=norm(xtemp(:)-W(:))/norm(err0_1);
    fprintf('Time elapsed is %f, residual is %d, original residual is %d\n',t,errTol,errOri);
else
    fprintf('Time elapsed is %f, residual is %d\n',t,errTol);
end

% figureObject(reshape(x0,m(1),m(2),NumElement),Z,m,NumElement,MU_e,1);
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
