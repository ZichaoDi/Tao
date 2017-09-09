x=linspace(min(a(:,2)),max(a(:,2)),100);
y=interp1(a(:,2),a(:,3),x);
ab=[x' y'];
for i=2:size(ab,1)-1
    x=ab(i-1:i+1,1);y=ab(i-1:i+1,2);
    K(i-1) = 2*abs((x(2)-x(1)).*(y(3)-y(1))-(x(3)-x(1)).*(y(2)-y(1))) ./ ...
      sqrt(((x(2)-x(1)).^2+(y(2)-y(1)).^2)*((x(3)-x(1)).^2+(y(3)-y(1)).^2)*((x(3)-x(2)).^2+(y(3)-y(2)).^2));
end

[~,ind]=max(K);
d=sum((repmat(ab(ind+1,:),size(a,1),1)-a(:,2:3)).^2,2);
[~,ind1]=min(d);
figure, plot(a(:,2),a(:,3),'m.-',ab(:,1),ab(:,2),'b--',ab(ind+1,1),ab(ind+1,2),'ro',a(ind1,2),a(ind1,3),'g*');




