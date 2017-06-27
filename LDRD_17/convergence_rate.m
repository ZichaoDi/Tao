function [q,l]=convergence_rate(x)
n=size(x,2);
for i=1:n-3
    num=norm(x(:,i+3)-x(:,i+2))/norm(x(:,i+2)-x(:,i+1));
    den=norm(x(:,i+2)-x(:,i+1))/norm(x(:,i+1)-x(:,i));
    q(i)=log(abs(num))/log(abs(den));
    l(i)=abs(num)/abs(den);
end 
