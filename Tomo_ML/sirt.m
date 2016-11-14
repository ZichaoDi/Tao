function [x,sing_rec,CR]=sirt(A,x,b,n_it)
% Input: sparse system matrix A, data b.
% Output: SIRT reconstruction x.
global N WS Ntot
ws=WS(:);
[rows cols] = size(A);
C = sparse(1 : cols, 1 : cols, 1 ./ sum(A));
R = sparse(1 : rows, 1 : rows, 1 ./ sum(A'));
CATR = C * A' * R;
sing_rec=[];
err=[];
res=[];
if(cols==N(1)^2)
   sing_rec(1,:)=[0,norm(b-A*x),norm(ws(:)-x)];
else
   sing_rec(1,:)=[0,norm(b-A*x),norm(x-x)];
end
fprintf('Iter    Residual          Error\n');
fprintf('%d      %e      %e\n',sing_rec(1,:));
for ki = 1 : n_it
    if(cols==N(1)^2)
        xOld=ws;
    else
        xOld=x;
    end
    x = x + CATR * (b - A * x);
    res(ki)=norm(A*x-b);
    err(ki)=norm(x-xOld);
    sing_rec(ki+1,:)=[ki,res(ki),err(ki)];
    if(mod(ki,100)==0 || ki==n_it)
        fprintf('%d      %e      %e\n',sing_rec(ki+1,:));
    end
    if(res(ki)<=1e-6)
        break;
    end
end
CR=(err(end)/err(1))^(1/ki);



