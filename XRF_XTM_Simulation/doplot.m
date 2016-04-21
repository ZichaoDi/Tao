function doplot(it,v,W_level,xinitial)
global current_n NumElement Z Element
%--------------------------------------------------------------
% Plot results of soap-film minimal surface
%--------------------------------------------------------------
m = [current_n current_n];
ws=reshape(W_level{1},prod(m)*NumElement,1);
xstar=v;
if(length(W_level)==2)
    nrow=5;
    xstarJ=W_level{2};
else
    nrow=4;
    ncol = NumElement;
end
%--------------------------------------------------------------
% set up boundary conditions
%--------------------------------------------------------------
 figure;
    clims=[0 max(ws)];
    for i=1:ncol
        subplot(nrow,ncol,i);
        
        errCom=reshape(xinitial(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-ws(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(errCom);colormap jet
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',...
        [],'YTick',[]);
        % colorbar;
        if(i==1)
            ylabel('Initial Guess','fontsize',12)
        end
        title(Element{Z(i)},'fontsize',12);
    
        subplot(nrow,ncol,i+1*ncol);
        
        errCom=reshape(xstar(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-ws(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(errCom);colormap jet
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',...
[],'YTick',[]);
        % colorbar;
        if(i==1)
            ylabel('XRF','fontsize',12)
        end

    
        subplot(nrow,ncol,i+2*ncol);
        
        errCom=reshape(ws(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));
        imagesc(errCom);colormap jet
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',...
        [],'YTick',[]);
        % colorbar;
        if(i==1)
            ylabel('Joint','fontsize',12)
        end
%         if(i==NumElement)
%         hp4 = get(subplot(nrow,NumElement,i+2*NumElement),'Position');
%         colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(3)]);
%         end
         subplot(nrow,ncol,i+3*ncol);
         ind=prod(m)*i-prod(m)+1:prod(m)*i;
         [~,sortInd]=sort(ws(ind));
         plot(1:prod(m),xinitial(sortInd),'ro',1:prod(m),xstar(sortInd),'bs',1:prod(m),ws(sortInd),'g*')
         xlim([0 prod(m)]);
         if(i==1)
             hleg=legend('initial','final','optimal','FontSize',6, 'Box', 'off','Orientation','horizontal');
             set(hleg,'units','pixels');
             lp=get(hleg,'outerposition');
             set(hleg,'Location','NorthWest', 'Box', 'off','outerposition',[lp(1)+50,15,lp(3),20]);
             ylabel('solution','fontsize',12)
         end
        if(nrow==5)
            subplot(nrow,ncol,i+4*ncol);
            
            errCom=reshape(xstarJ(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-ws(prod(m)*i-prod(m)+1:prod(m)*i
            imagesc(errCom);colormap jet
            set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',...
    [],'YTick',[]);
            % colorbar;
            if(i==1)
                ylabel('Joint','fontsize',12)
            end
        end
    end
%-------------------------------------------------------
T = sprintf('MG %1i',it);
% title(T)
