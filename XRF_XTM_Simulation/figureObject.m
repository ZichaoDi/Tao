function figureObject(W,Z,m,NumElement,MU_e,status)
persistent Map
O=zeros(m(1),m(2),NumElement); 
for e=1:NumElement
    O(:,:,e)=W(:,:,e).*Z(e);
end
Otest=sum(W.*repmat(reshape(Z,1,1,NumElement),m(1),m(2)),3);
MU=zeros(m);
for i=1:m(1)
    for j=1:m(2)
MU(i,j)=sum(reshape(W(i,j,:),NumElement,1).*reshape(MU_e(:,1,1),NumElement,1));
    end
end
if(status==0)
Map=figure('name','Attenuation Map');
% elseif(status==1)
%     figure('name','Initial Guess');
% elseif(status==2)
%     figure('name','Optimized Object')
end
% for e=1:NumElement
%     subplot(1,NumElement+1,e);
%     drawnow;
%     image(O(:,:,e));
%     axis xy
% end
% subplot(1,NumElement+1,NumElement+1);
figure(Map);
clims = [40 80];
subplot(1,5,status+1);
% imagesc(Otest,clims);colormap(gray);axis xy
if(status==0)
    title('Original');
elseif(status==1)
    title('Initial');
elseif(status==2);
    title('Joint');
    elseif(status==3);
    title('Diff');
    elseif(status==4);
    title('Alternate');
end
% subplot(1,2,2);
% drawnow;
% imagesc(MU); axis xy
% title('Attenuation Map at Incident Energy')