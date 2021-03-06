
omegaLocal=[-N/2 N/2 -N/2 N/2];
% x=linspace(omegaLocal(1),omegaLocal(2),N);
% y=linspace(omegaLocal(3),omegaLocal(4),N);
% [X,Y] = meshgrid((x(1:end-1)+x(2:end))./2,(y(1:end-1)+y(2:end))./2);
% I=(X-center(1)).^2+(Y-center(2)).^2;
C1=[N/4 N/4];
C2=[-N/4 -N/4];
C3=[N/4 -N/4];
r1=N/5; r2=N/6;
S1=[-N/4-r1 -N/4+r1 N/4-r1 N/4+r2];
S2=[C3(1)-r2,C3(1)+r2,C3(2)-r2, C3(2)+r2];
poly1=[C2(1) C2(2)-r1; C2(1)+r1 C2(2);C2(1) C2(2)+r1;C2(1)-r1 C2(2)];

mLocal=[N N]-1;
ne_local=5;
W=zeros(mLocal(1)+1,mLocal(2)+1,ne_local);
X = reshape(getNodalGrid(omegaLocal,mLocal),prod(mLocal+1),2);
p1=find(inCircle(X,repmat([C1,r1],prod(mLocal+1),1)));
p2=find(inCircle(X,repmat([C2,r1],prod(mLocal+1),1)));
p3=find(inCircle(X,repmat([C3,r2],prod(mLocal+1),1)));
p4=find(X(:,1)>=S1(1) & X(:,1)<=S1(2) & X(:,2)>=S1(3) & X(:,2)<=S1(4));
p5=find(X(:,1)>=S2(1) & X(:,1)<=S2(2) & X(:,2)>=S2(3) & X(:,2)<=S2(4));
p6=find(isPointInPolygon(X, poly1));
p7=setdiff(1:prod(mLocal+1),[p1;p2;p3;p4;p5]);
[dx,dy]=ind2sub(mLocal+1,[p7';p6]);
W(sub2ind([mLocal(1)+1,mLocal(2)+1,ne_local],dx,dy,1*ones(size(dx))))=0.6;
W(sub2ind([mLocal(1)+1,mLocal(2)+1,ne_local],dx,dy,2*ones(size(dx))))=1.2;
W(sub2ind([mLocal(1)+1,mLocal(2)+1,ne_local],dx,dy,3*ones(size(dx))))=0.2;

[indx,indy]=ind2sub(mLocal+1,setdiff([p4;p1;p2;p5],[p6;p3]));
W(sub2ind([mLocal(1)+1,mLocal(2)+1,ne_local],indx,indy,2*ones(size(indx))))=1.05;
W(sub2ind([mLocal(1)+1,mLocal(2)+1,ne_local],indx,indy,3*ones(size(indx))))=0.35;
W(sub2ind([mLocal(1)+1,mLocal(2)+1,ne_local],indx,indy,4*ones(size(indx))))=1.4;
W(sub2ind([mLocal(1)+1,mLocal(2)+1,ne_local],indx,indy,5*ones(size(indx))))=0.7;
[indx1,indy1]=ind2sub(mLocal+1,p3);
W(sub2ind([mLocal(1)+1,mLocal(2)+1,ne_local],indx1,indy1,2*ones(size(indx1))))=1.05;
W(sub2ind([mLocal(1)+1,mLocal(2)+1,ne_local],indx1,indy1,3*ones(size(indx1))))=0.35;
W(sub2ind([mLocal(1)+1,mLocal(2)+1,ne_local],indx1,indy1,4*ones(size(indx1))))=1.05;
W(sub2ind([mLocal(1)+1,mLocal(2)+1,ne_local],indx1,indy1,5*ones(size(indx1))))=1.05;
W=sum(W,3);
plotElement=0;
if(plotElement)
    clims=[0 max(W(:))];
    EleOri=[19 31 26 46 50];
    figure('name','Sample');
    for i=1:5
        ax(i)=subplot(1,5,i);
        imagesc(W(:,:,i),clims);
        if(i==1)
            xlabel('\leftarrow  0.1mm \rightarrow','FontSize',13,'FontWeight','bold')
        end
        title(Element{Z(i)},'FontSize',18,'FontWeight','bold');
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
        axis image xy
        colormap(jet)
    end
    h=colorbar('SouthOutside');
    set(h, 'Position', [.125 .35 .50 .03]);
%     set(get(h,'title'),'String','Density (g/cm^{3})');
% ylabel(h,'Density (g/cm^{3})','FontSize',13,'FontWeight','bold')
    for i=1:5
        pos=get(ax(i), 'Position');
        set(ax(i), 'Position', [pos(1) 0.1+pos(2) pos(3) 0.8*pos(4)]);
    end;
    drawnow;
end







