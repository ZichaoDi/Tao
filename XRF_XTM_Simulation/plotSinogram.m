figure,
for i=1:3, 
    subplot(1,3,i);imagesc(squeeze(d(i,:,:)));
    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]); 
    if(i==1)                                                            
        xlabel('\tau','FontSize',30,'FontWeight','bold');
         ylabel('\theta','FontSize',30,'FontWeight','bold'); 
    end;
end  
