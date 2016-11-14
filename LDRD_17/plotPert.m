Ntest=2;
DisR=[];
L_cr=[];
for ind_cr=1:Ntest;do_setup;L_cr(:,:,ind_cr)=L;DisR(:,:,ind_cr)=DisR_Simulated;end
cr_trans=cr./repmat(dz,Ntest,1)+N/2;
% save('testCircle_smaller.mat','cr_trans','DisR','W','thetan');       
% cmap=[0 I0];
% figure, 
% for i=1:Ntest, 
%     subplot(Ntest,2,2*i-1);imagesc(W);axis xy image;
%     hold on; plot(cr_trans(i,1),cr_trans(i,2),'r.','Markersize',10);
%     set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
%     subplot(Ntest,2,2*i);imagesc(map1D(DisR(:,:,i)',cmap));axis xy image;
% end

