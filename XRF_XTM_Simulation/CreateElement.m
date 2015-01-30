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
if(length(val)<4)
    subind=[1 2];
end
val=val(i2(subind));
NumElement=3;
W=zeros(m(1),m(2),NumElement);
for i=1:length(subind)
    Ws=zeros(m(1),m(2));
    Ws(abs(A(:)-val(i))<tol)=val(i)+ExtraVal(i);
    W(:,:,i)=Ws+0.1;
end
if(plotElement)
    figure('name','Element Map')
    for i=1:NumElement
        subplot(1,NumElement,i);
        clims=[0 1];
        imagesc(W(:,:,i),clims);
    end
end