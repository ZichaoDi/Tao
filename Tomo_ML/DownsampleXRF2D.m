%% Downsample the XRF spectrum
slice=30;
numChannel=2000;
nTau=1750;
% thetan=zeros(73,1);
DisXRF=zeros(1,1,2000);%zeros(length(thetan),nTau+1,numChannel);
DisR=zeros(length(thetan),nTau+1);
mStep=1750;%size(XRF);
stepx=ceil(mStep(1)/(nTau+1));
for t=1:length(thetan)
load(['2xfm_0',num2str(135+t),'.mat']);
% thetan(t)=angle;
% for i=1:nTau+1
%         if(i==nTau+1);
%             v=mStep(1);
%         else
%             v=i*stepx;
%         end
%         DisXRF(t,i,:)=sum(sum(XRF(((i-1)*stepx+1):v ,slice,:),1),2);
%         DisR(t,i)=sum(XTM( ((i-1)*stepx+1):v ,slice));
DisXRF=DisXRF+sum(sum(XRF,1),2);
end
% end
DetChannel=DetChannel(1:numChannel);
% save(['data/2xfm0415/slice',num2str(slice),'_',num2str(50),'.mat'],'DetChannel','DisR','DisXRF','thetan','nTau');
%%%%%%%%%%%%%%%%%%%%%%============================================================================
% ProJ=[173:245];
% slice=20;
% load thetan
% % [~,in]=unique(thetan);
% % ProJ=ProJ(in);
% load data/2xfm1211_14_4/DisR;
% m=size(DisR);
% I=zeros(m(2),m(2),m(1));
% % load DisR
% % DisR=zeros(51,1750,length(thetan));
% % DisRs=zeros(1750,50,length(thetan));
% for slice=1:m(1)
% % for i=1:length(thetan)
% %         load(['data/2xfm1211_14_4/2xfm_',num2str(ProJ(i)),'.mat']); 
% %         DisR(slice,:,i)=XTM(:,slice);
% % end
% I(:,:,slice)=iradon(squeeze(DisR(slice,1:end-2,:)),thetan,'spline','none','default',m(2));
% imagesc(I(:,:,slice));
% pause(1);
% end
% A=reshape(XRF_d,numxbins*numybins,mStep(3));

% for i=1:nTau+1
%         if(i==nTau+1);
%             v=1750;
%         else
%             v=i*10;
%         end
%         DisXRF(t,i,:)=sum(sum(XRF( ((i-1)*stepx+1):v ,slice,: ),1),2);
%         
% end
% {
% numx=5; numy=2;
% count=0;
% for i=1:floor(numxbins/numx):numxbins
%     for j=1:numy
% count=count+1;
% subplot(numx,numy,count)
% v=reshape(XRF_d(i+6,j,:),2000,1);
% semilogy(v,'LineWidth',1);
% counts(count)=sum(v);
% axis tight
% xlim([101,1100])
% set(gca,'xticklabel',[],'yticklabel',[])
% end
% end

%}
