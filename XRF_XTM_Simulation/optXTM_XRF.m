
global low up penalty
global W0 current_n
global SigMa_XTM SigMa_XRF
global NumElement err0 fiter nit maxiter xinitial

%%%----------------------------Initialize dependent variables
do_setup;
level=1;
current_n = N(level);
W= W_level{level};
XRF=xrf_level{1};
DisR=xtm_level{1};
%%---- smooth data
% DisR=BlurGaussian(DisR);
% for i_theta=1:size(DisR,2)
% DisR(:,i_theta)=smooth(DisR(:,i_theta));
% end
%%--------------------------------
L=L_level{level};
GlobalInd=GI_level{level};
SelfInd=SI_level{level};
m=m_level(level,:);
nTau=nTau_level(1);
SigMa_XTM=SigmaT{level};
SigMa_XRF=SigmaR{level};
if(NoSelfAbsorption)
    fprintf(1,'====== No Self Absorption, Transmission Detector Resolution is %d\n',nTau);
else
    fprintf(1,'====== With Self Absorption, Transmission Detector Resolution is %d\n',nTau);
end
%%%----------------------------------------------------------------------
%%%============== Rescale MU_e to make unity contribution
penalty=0;
xs_ele=[];
for ele=1;%:length(slice)
  
if(Joint==-1)
     fctn=@(W)sfun_XTM(W,DisR,MU_e,I0,L,thetan,m,nTau,NumElement); %% on W
     %fctn=@(MU)sfun_XTM_tensor(DisR,MU,I0,L,m,nTau);% NumElement=1;%% on attenuation coefficients \miu
elseif(Joint==1)
    fctn=@(W)sfun_Tensor_Joint(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0);
%     fctn_f=@(W)func_Tensor_Joint(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0);
%     fctn_J=@(W)sfun_Tensor_Joint_Jacobian(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0);
else
    fctn=@(W)sfun_Tensor(W,XRF,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau);
%     fctn_par=@(W)sfun_Tensor_par(W,XRF,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau);
%     fctn1=@(W)sfun_AdiMat(W,XRF,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau);
end

%-----------------------------------------------------------------------
rng('default')
x0=WS(:); %+0.04*rand(size(WS(:)));%
%{
x0=[];
for ii=1:length(slice),
x0(:,:,ii) =flipud(flipud( iR(:,:,slice(ii)))');
end
x0=x0(:)+1e0*rand(size(WS(:)));
x0(find(x0<0))=0;
%x0=x0(:);%*1e-2+rand(size(x0(:)));%reshape(iR(:,:,slice),[size(iR,1)*size(iR,2),1]);%10^(-1)*rand(nTol,1)+WS(:);%ww(:);%;%
load xs18_300_13TN_full
x0=xstar;
%}
xinitial=x0;
err0=norm(x0-W0);
e=cputime;
low=0*ones(size(x0));
up=inf*ones(size(x0));
% %%========================================================================
% [xstar,f,g,ierror] = tnbcm (x0,fctn,low,up,maxiter);
maxiter=100;
bounds=0;
if(bounds)
 [xstar,f,g,ierror] = tnbc (x0,fctn,low,up);
  save(['xs_Joint',num2str(N(1)),'_',num2str(maxiter),'_',num2str(numThetan),'TNBC_glassrod.mat'],'x0','xstar');
else
%  [xstar,f,g,histout,costdata] = gaussn(x0,fctn,1e-18,maxiter);
%  save(['xs',num2str(N(1)),'_',num2str(maxiter),'_',num2str(numThetan),'GN.mat'],'x0','xstar');
 [xstar,f,g,ierror] = tn (x0,fctn);
 save(['xs_Joint_real',num2str(N(1)),'_',num2str(maxiter),'_',num2str(numThetan),'tn_glassrod.mat'],'x0','xstar');
end

xs_ele(:,ele)=xstar;
end
%%%====================================================== Report Result

% for i=1:NumElement
%     err(i)=norm(xstar(prod(m)*i-(prod(m)-1):prod(m)*i)-WS(prod(m)*i-(prod(m)-1):prod(m)*i))/norm(xinitial(prod(m)*i-(prod(m)-1):prod(m)*i)-WS(prod(m)*i-(prod(m)-1):prod(m)*i));
% end
% figure('name','Elemental Residule')
% semilogy(1:NumElement,err,'r.-');
doplot(1,xstar, W_level);
if(Joint==1)
convFac_J=(fiter(end)/fiter(1))^(1/(nit+1));
t_J=cputime-e;
errTol_J=norm(xstar-WS(:))/norm(err0);
elseif(Joint==0)
  convFac_XRF=(fiter(end)/fiter(1))^(1/(nit+1));
  t_XRF=cputime-e;
errTol_XRF=norm(xstar-WS(:))/norm(err0);
end
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
     clims=[0 max(WS(:))];
    for i=1:NumElement
        subplot(4,NumElement,i);
        
        errCom=reshape(xinitial(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-WS(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(errCom,clims);colormap jet
        if(i==1)
            ylabel('Initial Guess','fontsize',12)
        end
        title(Element{Z(i)},'fontsize',12);
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
    end
    
    for i=1:NumElement
        subplot(4,NumElement,i+1*NumElement);
        
        errCom=reshape(xstar(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-WS(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(reshape(xstar,N,N));colormap jet
        if(i==1)
            ylabel('Final Soluction','fontsize',12)
        end
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
    end
    
    for i=1:NumElement
        subplot(4,NumElement,i+2*NumElement);
        
        errCom=reshape(WS(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));
        imagesc(errCom);colormap jet
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
        plot(1:prod(m),sort(xinitial(prod(m)*i-prod(m)+1:prod(m)*i)),'ro',1:prod(m),sort(xstar(prod(m)*i-prod(m)+1:prod(m)*i)),'bs',1:prod(m),sort(WS(prod(m)*i-prod(m)+1:prod(m)*i)),'g*')
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
