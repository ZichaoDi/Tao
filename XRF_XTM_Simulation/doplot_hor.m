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
    nrow=2;
    ncol = 3;
end
%--------------------------------------------------------------
% set up boundary conditions
%--------------------------------------------------------------
 figure;
    clims=[0 max(ws)];
    for i=1:nrow
        subplot(nrow,ncol,i+2*(i-1));
        
        errCom=reshape(xinitial(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-ws(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(errCom);colormap jet
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',...
        [],'YTick',[]);
        % colorbar;
        if(i==1)
            ylabel('Initial Guess','fontsize',12)
        end
        title(Element{Z(i)},'fontsize',12);
    
        subplot(nrow,ncol,i+2*i);
        
        errCom=reshape(xstar(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-ws(prod(m)*i-prod(m)+1:prod(m)*i
        imagesc(errCom);colormap jet
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',...
[],'YTick',[]);
        if(i==1)
            ylabel('XRF','fontsize',12)
        end

    
        subplot(nrow,ncol,i+2*ncol);
        
        errCom=reshape(ws(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));
        imagesc(errCom);colormap jet
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',...
        [],'YTick',[]);
        if(i==1)
            ylabel('Joint','fontsize',12)
        end
    end
%-------------------------------------------------------
