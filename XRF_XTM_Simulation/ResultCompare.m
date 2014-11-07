close all
colMap={'hot','gray','bone','pink','copper'};
Ele={'Mo','Cu','Fe','Ca','In'};
uplim=max([max(W(:)),max(W2(:)),max(W3(:))]);%,max(W1(:))
for i=1:NumElement
subplot(3,NumElement,i)

imagesc(W(:,:,i));
colormap(colMap{i})
set(gca,'xtickLabel',[],'ytickLabel',[],'clim',[0 uplim],'FontSize',15,'FontWeight','bold');
axis square, box on;
title([Ele{i}])
set(gca,'FontSize',10,'FontWeight','bold');
if(i==1)
    ylabel('True')
end

freezeColors
end

% for i=1:NumElement
% subplot(4,NumElement,NumElement+i)
% imagesc(W1(:,:,i));
% colormap(colMap{i})
% set(gca,'xtickLabel',[],'ytickLabel',[],'clim',[0 uplim]);
% 
% axis square, box on;
% set(gca,'FontSize',10,'FontWeight','bold');
% if(i==1)
%     ylabel('JRT')
% end
% freezeColors
% end

for i=1:NumElement
subplot(3,NumElement,1*NumElement+i)
imagesc(W2(:,:,i));
colormap(colMap{i})
set(gca,'xtickLabel',[],'ytickLabel',[],'clim',[0 uplim]);
axis square, box on;
if(i==1)
    ylabel('SRT')
end
freezeColors
end

for i=1:NumElement
subplot(3,NumElement,2*NumElement+i)
imagesc(W3(:,:,i));
colormap(colMap{i})
set(gca,'xtickLabel',[],'ytickLabel',[],'clim',[0 uplim]);

axis square, box on;
if(i==1)
    ylabel('ART')
end
freezeColors
end
% saveas(gca,['/Users/Wendydi/Dropbox/TaoFusion/presentation/images/ResultCom50_partial.png']) 
print('-depsc','-tiff','-r300','./result/ResultCom50_partial');