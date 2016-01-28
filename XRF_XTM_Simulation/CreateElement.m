global NumElement

A=phantom(m(1));
% load Logan50;
% A=Logan;
% A=ones(m);
A(A(:)<0)=0;
tol=eps^(1/2);
val=unique(A(:));
% val=val(val~=0);
val_i=[];
for i_v=1:length(val)
    val_i(i_v)=length(find(A(:)==val(i_v)));
end
ExtraVal=[0.1 0.2 0.3];
[i1,i2]=sort(val_i,'descend');
subind=[1 2 4];
subind=subind(1:NumElement);
if(length(val)<4)
    subind=[1 2];
    if(NumElement==1)
        subind=1;
    end
end
val=val(i2(subind));
W=zeros(m(1),m(2),NumElement);
for i=1:length(subind)
    Ws=zeros(m(1),m(2));
    Ws(abs(A(:)-val(i))<tol)=val(i)+ExtraVal(i);
    W(:,:,i)=Ws;%+0.1;
end
if(plotElement)
    figure('name','Element Map')
    clims=[0 max(W(:))];
    for i=1:NumElement
        ax(i)=subplot(1,NumElement,i);
        imagesc(W(:,:,i),clims);
        if(i==1)
            xlabel('\leftarrow       0.1mm      \rightarrow','FontSize',13,'FontWeight','bold')
        end
        title(Element{Z(i)},'FontSize',18,'FontWeight','bold');
        set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
        axis image xy
        colormap(jet)
    end
    h=colorbar('SouthOutside');
    set(h, 'Position', [.125 .28 .50 .03]);
    for i=1:NumElement
        pos=get(ax(i), 'Position');
        set(ax(i), 'Position', [pos(1) 0.1+pos(2) pos(3) 0.8*pos(4)]);
    end;
    drawnow;
end
