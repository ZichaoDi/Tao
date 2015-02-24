
% global gv
global low up penalty gama
global W0 current_n N
global SigMa_XTM SigMa_XRF
global fctn_f err0

%%%----------------------------Initialize dependent variables
level=find(N==current_n);
W= W_level{level};
XRF=xrf_level{1};
DisR=xtm_level{1};
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
    fctn=@(W)sfun_Tensor_Joint(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0);
%     fctn_f=@(W)func_Tensor_Joint(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0);
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

fctn_J=@(W)sfun_Tensor_Joint_Jacobian(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0);
feval(fctn_J,x0);
% tic;
% [f,g]=feval(fctn,W(:));
% toc;
% return;
err0=norm(x0-W0);
ws=Wtest(:);
e=cputime;
low=0*ones(size(x0));
up=1e6*ones(size(x0));
% %%========================================================================
%  options = optimset('Algorithm','interior-point','Display','iter','MaxFunEvals',10000);%,'DerivativeCheck','off','Diagnostics','off','GradObj','off','Display','iter','AlwaysHonorConstraints','none','TolCon',1e-10,'TolX',1e-16,'TolFun',1e-15);%
%  [xstar,fval] = fmincon(fctn,x0,[],[],[],[],zeros(size(x0)),[],[],options);
%  return;
% [xstar,f,g,ierror] = tnbcm (x0,fctn,low,up,maxiter);
% options = optimset('Display','iter','TolFun',1e-8);
% xstar = lsqnonlin(fctn,x0,[],[],options);
[xstar,f,g,ierror] = tnbc (x0,fctn,low,up);
%   [xstar,f,g,ierror] = tn (x0,fctn);
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
% if(DiscreteScale)
%     AbsErr=norm(xtemp(:)-W(:))
%     IniErr=norm(err0_1)
%     errOri=norm(xtemp(:)-W(:))/norm(err0_1);
%     fprintf('Time elapsed is %f, residual is %d, original residual is %d\n',t,errTol,errOri);
% else
%     fprintf('Time elapsed is %f, residual is %d\n',t,errTol);
% end

% figureObject(reshape(x0,m(1),m(2),NumElement),Z,m,NumElement,MU_e,1);
if(plotResult)
    figure('name','Elemental Residule');
    %     clims=[0 max([xinitial;ws;xstar])];
    for i=1:NumElement
        subplot(4,NumElement,i);
        
        errCom=reshape(xinitial(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-ws(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(errCom,clims);colormap jet
        if(i==1)
            ylabel('Initial Guess','fontsize',12)
        end
        title(Element{Z(i)},'fontsize',12);
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
    end
    
    for i=1:NumElement
        subplot(4,NumElement,i+1*NumElement);
        
        errCom=reshape(xstar(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-ws(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(errCom,clims);colormap jet
        if(i==1)
            ylabel('Final Soluction','fontsize',12)
        end
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
    end
    
    for i=1:NumElement
        subplot(4,NumElement,i+2*NumElement);
        
        errCom=reshape(ws(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));
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
        plot(1:prod(m),sort(xinitial(prod(m)*i-prod(m)+1:prod(m)*i)),'ro',1:prod(m),sort(xstar(prod(m)*i-prod(m)+1:prod(m)*i)),'bs',1:prod(m),sort(ws(prod(m)*i-prod(m)+1:prod(m)*i)),'g*')
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
% %%%%%%%%%%%%%%%%%%%%=======================================================
%%%===================================================== Derivative Test
%   gh=foo(fctn,x0);
% % % % % foo(fctn1,x0);
% tic;
%  [f,g]=feval(fctn,W(:));
% T_Tensor=toc;
% tic;
%  [f1,g1]=feval(fctn1,W(:));
% T_AD=toc;
% return;
% gh=foo(fctn,x0);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% For loop version (previous)
% % disp('========= For Loop Version')
% % TolP=3;
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
