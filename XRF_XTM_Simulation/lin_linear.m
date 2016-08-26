function [x,alpha]=lin_linear(e,x,A,b);
alpha=1;
for i=1:10
    x_new=x+alpha*e;
    if(norm(A*x_new-b)<norm(A*x-b))
        x=x_new;
        break;
    else
        alpha=alpha/2;
    end
end
fprintf('step length is %d\n',alpha);
