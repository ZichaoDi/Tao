function [f,g]=Quadratic(x,L)

f=1/2*x'*L*x-sum(x);
g=L*x-1;
