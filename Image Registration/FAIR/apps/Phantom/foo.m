
function foo(fctn,x0)
z0 = x0;
z = z0(:);
n = length(z);
h = (eps)^(1/2);
[f,para,g] =feval(fctn,z0);
% [f,g] =feval(fctn,z0);
gh = g;

diff_ind = 1; % Use 0 for central differencing, 1 for forward differencing

if (diff_ind == 1);
    disp('Forward Differencing')
    fprintf(' i                g                gh         error\n')
    for i= 1:n;
        zh = z;
        zh(i) = zh(i) + h;
         [fh,para,gh] = feval(fctn,zh);
%  [fh,gh] = feval(fctn,zh);
        gh(i) = (fh-f)/h;
        eee = abs(g(i)-gh(i))/(1+abs(g(i)));
        fprintf('%3i     %18.8e        %18.8e     %8.1e\n',i,g(i),gh(i),eee);
    end;
    err = norm(g-gh)/(1+norm(g))
else
    disp('Central Differencing')
    fprintf(' i               fh1                  fh2                  g                       gh              error\n')
    for i= 1:n;
        zh = z;
        zh(i) = zh(i) + h;
        fh1 = feval(fctn,zh);
        zh = z;
        zh(i) = zh(i) - h;
        fh2 = feval(fctn,zh);
        gh(i) = (fh1-fh2)/(2*h);
        eee = abs(g(i)-gh(i))/(1+abs(g(i)));
        fprintf('%3i     %18.8e  %18.8e  %18.8e  %18.8e     %8.1e\n',i,fh1,fh2,g(i),gh(i),eee);
    end;
    err = norm(g-gh,'inf')/(1+norm(g,'inf'))
end;