function [f,g]=linear(x,A,b);

f=sum((A*x-b).^2);
g=2*A'*A*x-2*A'*b;
