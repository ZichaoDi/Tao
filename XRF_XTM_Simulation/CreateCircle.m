close all;
Dis=800;
omega=[-Dis/2 Dis/2 -Dis/2 Dis/2];
% x=linspace(omega(1),omega(2),Dis);
% y=linspace(omega(3),omega(4),Dis);
% [X,Y] = meshgrid((x(1:end-1)+x(2:end))./2,(y(1:end-1)+y(2:end))./2);
% I=(X-center(1)).^2+(Y-center(2)).^2;
C1=[Dis/4 Dis/4];
C2=[-Dis/4 -Dis/4];
C3=[Dis/4 -Dis/4];
r1=Dis/5; r2=Dis/6;
S1=[-Dis/4-r1 -Dis/4+r1 Dis/4-r1 Dis/4+r2];
S2=[C3(1)-r2,C3(1)+r2,C3(2)-r2, C3(2)+r2];
poly1=[C2(1) C2(2)-r1; C2(1)+r1 C2(2);C2(1) C2(2)+r1;C2(1)-r1 C2(2)];

m=[Dis Dis];
NumElement=5;
W=zeros(m(1)+1,m(2)+1,NumElement);
X = reshape(getNodalGrid(omega,m),prod(m+1),2);
p1=find(inCircle(X,repmat([C1,r1],prod(m+1),1)));
p2=find(inCircle(X,repmat([C2,r1],prod(m+1),1)));
p3=find(inCircle(X,repmat([C3,r2],prod(m+1),1)));
p4=find(X(:,1)>=S1(1) & X(:,1)<=S1(2) & X(:,2)>=S1(3) & X(:,2)<=S1(4));
p5=find(X(:,1)>=S2(1) & X(:,1)<=S2(2) & X(:,2)>=S2(3) & X(:,2)<=S2(4));
p6=find(isPointInPolygon(X, poly1));
size(p6)
p7=setdiff(1:prod(m+1),[p1;p2;p3;p4;p5]);
[dx,dy]=ind2sub(m+1,[p7';p6]);
W(sub2ind([m(1)+1,m(2)+1,NumElement],dx,dy,1*ones(size(dx))))=0.6;
W(sub2ind([m(1)+1,m(2)+1,NumElement],dx,dy,2*ones(size(dx))))=1.2;
W(sub2ind([m(1)+1,m(2)+1,NumElement],dx,dy,3*ones(size(dx))))=0.2;

[indx,indy]=ind2sub(m+1,setdiff([p4;p1;p2;p5],[p6;p3]));
W(sub2ind([m(1)+1,m(2)+1,NumElement],indx,indy,2*ones(size(indx))))=1.05;
W(sub2ind([m(1)+1,m(2)+1,NumElement],indx,indy,3*ones(size(indx))))=0.35;
W(sub2ind([m(1)+1,m(2)+1,NumElement],indx,indy,4*ones(size(indx))))=1.4;
W(sub2ind([m(1)+1,m(2)+1,NumElement],indx,indy,5*ones(size(indx))))=0.7;
[indx1,indy1]=ind2sub(m+1,p3);
W(sub2ind([m(1)+1,m(2)+1,NumElement],indx1,indy1,2*ones(size(indx1))))=1.05;
W(sub2ind([m(1)+1,m(2)+1,NumElement],indx1,indy1,3*ones(size(indx1))))=0.35;
W(sub2ind([m(1)+1,m(2)+1,NumElement],indx1,indy1,4*ones(size(indx1))))=1.05;
W(sub2ind([m(1)+1,m(2)+1,NumElement],indx1,indy1,5*ones(size(indx1))))=1.05;

figure('name','Sample');
for i=1:5
subplot(2,3,i),
clims=[0 2];
imagesc(W(:,:,i),clims);
axis image xy
colormap(jet)
end









