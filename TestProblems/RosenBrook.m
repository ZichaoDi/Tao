function [f,g]=RosenBrook(x)
n=length(x);
f=sum(100*(x(2:end)-x(1:end-1).^2).^2+(x(1:end-1)-1).^2);
g=zeros(size(x));
for i=1:n
    if(i==1)
        g(i)=-400*(x(i+1)-x(i)^2)*x(i)+2*(x(i)-1);
    elseif(i==n)
        g(i)=200*(x(i)-x(i-1)^2);
    else
        g(i)=200*(x(i)-x(i-1)^2)-400*(x(i+1)-x(i)^2)*x(i)+2*(x(i)-1);
    end
end