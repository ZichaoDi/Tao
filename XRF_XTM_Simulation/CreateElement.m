
% A=phantom(m(1));
% load Logan50;
% A=Logan;
A=ones(m);
A(A(:)<0)=0;
tol=eps^(1/2);
val=unique(A(:));
val=val(val~=0);
NumElement=length(val);
W=zeros(m(1),m(2),NumElement);

for i=1:NumElement
Ws=zeros(m(1),m(2));
Ws(abs(A(:)-val(i))<tol)=val(i)+1;
W(:,:,i)=Ws+0.1;
end
if(plotElement)
 figure('name','Element Map')
 for i=1:NumElement
subplot(1,NumElement,i);imagesc(W(:,:,i));
 end
end