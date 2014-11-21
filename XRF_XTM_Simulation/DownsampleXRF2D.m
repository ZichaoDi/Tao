%% Downsample the XRF spectrum
slice=20;
numChannel=2000;
nTau=50-1;
thetan=zeros(36,1);
DisXRF=zeros(length(thetan),nTau+1,numChannel);
DisR=zeros(length(thetan),nTau+1);
mStep=size(XRF);
stepx=ceil(mStep(1)/(nTau+1));
for t=1:length(thetan)
load(['2xfm_0',num2str(135+t),'.mat']);
thetan(t)=angle;
for i=1:nTau+1
        if(i==nTau+1);
            v=mStep(1);
        else
            v=i*stepx;
        end
        DisXRF(t,i,:)=sum(sum(XRF(((i-1)*stepx+1):v ,slice,:),1),2);
        DisR(t,i)=sum(XTM( ((i-1)*stepx+1):v ,slice));
end
end
DetChannel=DetChannel(1:numChannel);
save(['data/2xfm0415/slice',num2str(slice),'_',num2str(50),'.mat'],'DetChannel','DisR','DisXRF','thetan','nTau');
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
%{
numx=5; numy=2;
count=0;
for i=1:floor(numxbins/numx):numxbins
    for j=1:numy
count=count+1;
subplot(numx,numy,count)
v=reshape(XRF_d(i+6,j,:),2000,1);
semilogy(v,'LineWidth',1);
counts(count)=sum(v);
axis tight
xlim([101,1100])
set(gca,'xticklabel',[],'yticklabel',[])
end
end

%}