load ComA
A=Acrop;
FS=16;
close all
nm=1;
nt=1;%linspace(0,4,nm);
ee=[];
for i=1:nm
n=nt(i);
 x1=log(A(:,1).^n).^n; x2=log(A(:,2));y=log(A(:,3));
% x1=A(:,1).^n; x2=A(:,2);y=A(:,3);

X = [ones(size(x1)) x1 x2 x1.*x2];
[b,bint,r] = regress(y,X); 
ee(i)=norm(r);
figure,
scatter3(A(:,1),A(:,2),A(:,3),'filled')
hold on
x1fit = linspace(min(A(:,1)),max(A(:,1)),20);
x2fit = linspace(min(A(:,2)),max(A(:,2)),20);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = exp(b(1) + b(2)*log(X1FIT.^n).^n + b(3)*log(X2FIT) + b(4)*log(X1FIT.^n).^n.*log(X2FIT));
mesh(X1FIT,X2FIT,YFIT,'LineWidth',0.03)

% %%%%------------------------------------------------------
% figure,
% scatter3(x1,x2,y,'filled')
% hold on
% x1fit = linspace(min(x1),max(x1),20);
% x2fit = linspace(min(x2),max(x2),20);
% [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
% YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
% mesh(X1FIT,X2FIT,YFIT,'LineWidth',0.03)
end
ylabel('$|\Theta|$','FontSize',FS,'FontWeight','bold','Interpreter','latex')
xlabel('$|\mathcal{V}|$','FontSize',FS,'FontWeight','bold','Interpreter','latex')
zlabel('Time Elapsed (sec)','FontSize',FS,'FontWeight','bold')
set(gca,'FontSize',FS-2)
% [e,ei]=min(ee);
%%%%%%%%%%%%%%%%%%%%-------------------------------------------
theta=1:10;
v=unique(A(:,1)./3);
n=unique(A(:,1));
nn=length(n);
rng('default');
figure,
f=18.9*1e-4*repmat(v.^1.5,1,size(theta,2)).*repmat(theta,size(v,1),1);
cmap = colormap(rand(100,3));
nn_sub=[1:4,6:nn];
for ii=1:nn-1;%size(f,1)
    i=nn_sub(ii);
    semilogy(theta,f(i,:),'--','linewidth',2,'color',cmap(i+5,:));
        hold on;
    ind=find(A(:,1)==n(i));
    semilogy(A(ind,2),A(ind,3),'o','linewidth',2,'color',cmap(i+5,:));
    text(A(ind(end),2)+0.2,A(ind(end),3),['$|\mathcal{V}|$=',num2str(v(i))],'FontSize',16,'Interpreter','latex','FontWeight','bold')
end

axis([1 12 min(A(:,3)) 1e3])
xlabel('$|\Theta|$','FontSize',FS,'FontWeight','bold','Interpreter','latex')
ylabel('Time Elapsed (sec)','FontSize',FS,'FontWeight','bold')
set(gca,'FontSize',FS-2)





% for i=2:10
% v=(A(:,2)==i);
% coeff=[ones(sum(v),1) log(A(v,1))]\log(A(v,3));
% end
% 
% figure
% v=(A(:,2)>0);
% AA=[ones(sum(v),1)  log(A(v,1))];
% coeff=AA\(log(A(v,3))-log(A(v,2)))
% plot(log(A(v,3))-log(A(v,2)),AA*coeff,'x') 


