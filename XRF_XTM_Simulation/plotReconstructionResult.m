close all
N=20;
load Phantom20;
w=Phantom20;
load PeriodicTable
xs=[];
load xs20X_10_noise
xsR=xs20X_10_noise;
load xs20J_10_noise
xsJ=xs20J_10_noise;
%%%-----------------------------------------
% load xs20_10X
% xsR=xs20_10X;%eval('xs10_6X.mat');
% load xs20_10J
% xsJ=xs20_10J;%eval('xs10_6J.mat');

% load xs20_10XRT
% xsT=xstar;

%%%-----------------------------------------
% load Phantom3;
% w=Phantom3;
% load xs3_4J
% load xs3_4X;%pNoSelf
% load xs3_4T;
Z=[19 31 26];% 46 50
NumElement=length(Z);
clims=[0,max(w(:))];
numRow=7;
%-------------------------------- plot results from different beta

% for i=1:NumElement+1,
%     if(i~=NumElement+1)
%         subplot(numRow,NumElement+1,i),imagesc(w(:,:,i),clims);
%         if(i==1)
%             ylabel('True Solution','FontSize',16,'FontWeight','bold')
%         end
%         title(Element{Z(i)},'FontSize',18,'FontWeight','bold');
%     else
%         subplot(numRow,NumElement+1,NumElement+1),imagesc(zeros(size(w,1),size(w,2)),clims);
%         title('Error','FontSize',18,'FontWeight','bold');
%     end
%     colormap(jet),axis image xy,
%     set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
%     
% end
% for j=1:6
%     load(['xsBeta',num2str(j),'.mat']);
%     xs=reshape(xstar,10,10,NumElement);
%     e=sum(abs(xs-w),3);
%     err=norm(e);
%     for i=1:NumElement+1,
%         if(i~=NumElement+1)
%             subplot(numRow,NumElement+1,j*(NumElement+1)+i),imagesc(xs(:,:,i),clims);
%             if(i==1)
%                 ylabel(['Solution',num2str(j)],'FontSize',16,'FontWeight','bold')
%             end
%             title(Element{Z(i)},'FontSize',18,'FontWeight','bold');
%         else
%             subplot(numRow,NumElement+1,j*(NumElement+1)+i),imagesc(e,clims);
%             title('Error','FontSize',18,'FontWeight','bold');
%         end
%         colormap(jet),axis image xy,
%         set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
%         
%     end
% end
%-------------------------------------------------------------------
rng('default')

w0=w+0.1*rand(size(w));
 e0=sum(abs(w0-w),3);
w2=reshape(xsJ,N,N,NumElement);%w+0.06*rand(size(w));
 e2=sum(abs(w2-w),3);
w1=reshape(xsR,N,N,NumElement);%w+0.01*rand(size(w));
e1=sum(abs(w1-w),3);
% w3=reshape(xsT,N,N,NumElement);%w+0.01*rand(size(w));
% e3=sum(abs(w3-w),N);
NumRow=4;
for i=1:NumElement+1,
    if(i~=NumElement+1)
        subplot(NumRow,NumElement+1,i),imagesc(w(:,:,i),clims);
        if(i==1)
            ylabel('True Solution','FontSize',16,'FontWeight','bold')
        end
        title(Element{Z(i)},'FontSize',18,'FontWeight','bold');
    else
        subplot(NumRow,NumElement+1,NumElement+1),imagesc(zeros(size(w,1),size(w,2)),clims);
        title('Error','FontSize',18,'FontWeight','bold');
    end
    colormap(jet),axis image xy,
    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    
end
for i=1:NumElement+1,
    if(i~=NumElement+1)
        subplot(NumRow,NumElement+1,NumElement+1+i),imagesc(w0(:,:,i),clims);
        if(i==1)
            ylabel('Initial Guess','FontSize',16,'FontWeight','bold')
        end
    else
        subplot(NumRow,NumElement+1,2*(NumElement+1)),imagesc(e0,clims);
    end
    colormap(jet),axis image xy,
    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
end
for i=1:NumElement+1,
    if(i~=NumElement+1)
        
        subplot(NumRow,NumElement+1,2*(NumElement+1)+i),imagesc(w1(:,:,i),clims);
        if(i==1)
            ylabel('XRF','FontSize',16,'FontWeight','bold')
        end
    else
        subplot(NumRow,NumElement+1,3*(NumElement+1)),imagesc(e1,clims);
    end
    colormap(jet),axis image xy,
    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
end
for i=1:NumElement+1,
    if(i~=NumElement+1)
        
        subplot(NumRow,NumElement+1,3*(NumElement+1)+i),imagesc(w2(:,:,i),clims);
        if(i==1)
            ylabel('JRT','FontSize',16,'FontWeight','bold')
        end
    else
        subplot(NumRow,NumElement+1,4*(NumElement+1)),imagesc(e2,clims);
    end
    colormap(jet),axis image xy,
    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    
%     if(i~=NumElement+1)
%         
%         subplot(5,NumElement+1,4*(NumElement+1)+i),imagesc(w3(:,:,i),clims);
%         if(i==1)
%             ylabel('XRT','FontSize',11,'FontWeight','bold')
%         end
%     else
%         subplot(5,NumElement+1,5*(NumElement+1)),imagesc(e3,clims);
%     end
%     colormap(jet),axis image xy,
%     set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    
end
    h=colorbar('SouthOutside');
    set(h, 'Position', [.2 .05 .50 .02]);