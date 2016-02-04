function [f,g]=linear(x,A,b);

% f=sum((A*x-b).^2);
f=(A*x-b)'*(A*x-b);
g=2*A'*A*x-2*A'*b;
