function doplot(it,v,W_level)
global xinitial current_n NumElement Z Element WS
%--------------------------------------------------------------
% Plot results of soap-film minimal surface
%--------------------------------------------------------------
m = [current_n current_n];
%reshape(W_level{1},prod(m)*NumElement,1);
NumElement=1;
xstar=v;
W0=WS(:);
%--------------------------------------------------------------
% set up boundary conditions
%--------------------------------------------------------------
 figure(24);
    clims=[0 max(W0(:))];
    for i=1:NumElement
        subplot(4,NumElement,i);
        
        errCom=reshape(xinitial(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-W0(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(errCom,clims);colormap jet
        if(i==1)
            ylabel('Initial Guess','fontsize',12)
        end
        title(Element{Z(i)},'fontsize',12);
    end
    
    for i=1:NumElement
        subplot(4,NumElement,i+1*NumElement);
        
        errCom=reshape(xstar(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-W0(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(errCom,clims);colormap jet
        if(i==1)
            ylabel('Final Soluction','fontsize',12)
        end
    end
    
    for i=1:NumElement
        subplot(4,NumElement,i+2*NumElement);
        
        errCom=reshape(W0(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));
        imagesc(errCom,clims);colormap jet
        if(i==1)
            ylabel('True Soluction','fontsize',12)
        end
%         if(i==NumElement)
%         hp4 = get(subplot(4,NumElement,i+2*NumElement),'Position');
%         colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(3)]);
%         end
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
    end
%-------------------------------------------------------
T = sprintf('MG %1i',it);
title(T)
