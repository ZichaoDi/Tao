% global f_xrf f_xtm
% do_setup;
% DecomposedElement=0;
% Be=[0 0.25 0.5 1 2 4 8];
% xs=zeros(length(Be),N*N*NumElement);
% perato=[];
% for i_b=1:length(Be)
%     Beta=Be(i_b);
%     opt;
%     perato=[perato;[Beta f_xrf f_xtm]];xs(i_b,:)=xstar;
% end
% save (['perato_',sample, num2str(N),'_',num2str(beta_d), '.mat'],'perato','xs')

load perator_Rod195_10.mat
figure, 
for j=1:size(perato,1); 
    for i=1:NumElement, 
        xstar=reshape(xs(j,:),N,N,NumElement);subplot(size(perato,1), NumElement, NumElement*(j-1)+i); imagesc(xstar(:,:,i));
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
        if(i==1)
            ylabel(['\beta = ',num2str(Be(j))]);
        end
        if(j==1)
            title(Element)
    end;
end
[~,ind]=sort(perato(:,2)); figure, plot(perato(ind,2),perato(ind,3),'r.-');
xlabel('XRF residual');ylabel('XRT residual');
for j=1:size(perato,1),
    text(perato(ind(j),2),perato(ind(j),3),num2str(perato(ind(j),1)));
end 

