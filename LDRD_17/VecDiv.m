function [Dx,Dy] = VecDiv(n)
Tik=diag([-1;ones(n-2,1);1])+diag(1/2*ones(n-1,1),1)+diag(-1/2*ones(n-1,1),-1);
Dx=kron(eye(n),Tik);
Dy=sparse(n^2,n^2);

for i=1:n^2

    if(i<=n)
        Dy(i,i)=-1; Dy(i,i+n)=1;
    elseif(i> n^2-n)
        Dy(i,i-n)=-1; Dy(i,i)=1;
    else
        Dy(i,i-n)=-1; Dy(i,i+n)=-1;
    end
end





