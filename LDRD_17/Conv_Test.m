%% Convergence Factor Test
%% Phantom Test (3 chemical elements) with 1 angle
% do_setup;
% fid = fopen('Convergence_Factor.txt','a');
% for i=1:length(N)
%     current_n=N(i);
%     Joint=1;
%     optXTM_XRF_Tensor;
%     Joint=0;
%     optXTM_XRF_Tensor;
%     ntol=current_n^2*NumElement;
%     fprintf(fid,'%d    %12.4e     %12.4e    %12.4e     %12.4e    %12.4e     %12.4e\n', ntol, convFac_J, t_J, errTol_J, convFac_XRF, t_XRF, errTol_XRF);
% end
%
% fclose(fid);

close all
% load Phase20;
% w=Phase20;
% load PeriodicTable
% load xs20_12Jp
% load xs20_12Xp;%pNoSelf
load Phantom20;
w=Phantom20;
load PeriodicTable
load xs20_10J
load xs20_10X;%pNoSelf
Z=[19 31 26];%  46 50];%
NumElement=length(Z);
rng('default')
w0=0.1*rand(size(w));
 e0=abs(sum(w0-w,3));
w2=reshape(xs20_10J,20,20,NumElement);%w+0.06*rand(size(w));
 e2=abs(sum(w2-w,3));
w1=reshape(xs20_10X,20,20,NumElement);%w+0.01*rand(size(w));
e1=abs(sum(w1-w,3));
clims=[0,max(w(:))];

for i=1:NumElement+1,
    if(i~=NumElement+1)
        subplot(4,NumElement+1,i),imagesc(w(:,:,i),clims);
        if(i==1)
            ylabel('True Solution','FontSize',11,'FontWeight','bold')
        end
        title(Element{Z(i)},'FontSize',18,'FontWeight','bold');
    else
        subplot(4,NumElement+1,NumElement+1),imagesc(zeros(size(w,1),size(w,2)),clims);
        title('Error','FontSize',18,'FontWeight','bold');
    end
    colormap(jet),axis image xy,
    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    
end
for i=1:NumElement+1,
    if(i~=NumElement+1)
        subplot(4,NumElement+1,NumElement+1+i),imagesc(w0(:,:,i),clims);
        if(i==1)
            ylabel('Initial Guess','FontSize',11,'FontWeight','bold')
        end
    else
        subplot(4,NumElement+1,2*(NumElement+1)),imagesc(e0,clims);
    end
    colormap(jet),axis image xy,
    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
end
for i=1:NumElement+1,
    if(i~=NumElement+1)
        
        subplot(4,NumElement+1,2*(NumElement+1)+i),imagesc(w1(:,:,i),clims);
        if(i==1)
            ylabel('AFT','FontSize',11,'FontWeight','bold')
        end
    else
        subplot(4,NumElement+1,3*(NumElement+1)),imagesc(e1,clims);
    end
    colormap(jet),axis image xy,
    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
end
for i=1:NumElement+1,
    if(i~=NumElement+1)
        
        subplot(4,NumElement+1,3*(NumElement+1)+i),imagesc(w2(:,:,i),clims);
        if(i==1)
            ylabel('JFT','FontSize',11,'FontWeight','bold')
        end
    else
        subplot(4,NumElement+1,4*(NumElement+1)),imagesc(e2,clims);
    end
    colormap(jet),axis image xy,
    set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
end
