%% Downsample the XRF spectrum
load 2xfm_0209.mat
mStep=size(XRF);
numxbins=19;
numybins=2;
stepy=ceil(mStep(2)/numybins);
stepx=ceil(mStep(1)/numxbins);
XRF_d=zeros(numxbins,numybins,mStep(3));
XTM_d=zeros(numxbins,numybins);

for i=1:numxbins
    for j=1:numybins
%         if j==1
%             v=1:26;
%         elseif j==2 
%             v=26:51;
%         end
if(j==numybins)
    v=mStep(2);
else
    v=j*stepy;
end
        XRF_d(i,j,:)=sum(sum(XRF( ((i-1)*stepx+1):(i*stepx) ,((j-1)*stepy+1):v ,: ),2),1);
        XTM_d(i,j)=sum(sum(XTM( ((i-1)*stepx+1):(i*stepx) ,((j-1)*stepy+1):v)));
    end
end



% A=reshape(XRF_d,numxbins*numybins,mStep(3));


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