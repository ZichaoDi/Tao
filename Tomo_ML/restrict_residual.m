function BB=restrict_residual(n)
global coarsen_type
if(strcmp(coarsen_type,'oscillate'))
    m=ceil(n/3);
    BB=sparse(zeros(m,n));
    for i=1:m
        if(i==m)
            if(3*i-2==n)
                BB(m,3*m-2)=1;
            else
                BB(m,3*m-2)=1;
                BB(m,3*m-1)=2;
            end
        else
                BB(i,3*i-2)=1;
                BB(i,3*i-1)=2;
                BB(i,3*i)=1;
        end
    end
elseif(strcmp(coarsen_type,'smooth'))
    if(mod(n,2)==0)
        m=n/2;
    else
        m=(n+1)/2;
    end
    BB=sparse(zeros(m,n));
    for i=1:m
        if(i==1)
                BB(i,2*i-1)=2;
                BB(i,2*i)=1;
        elseif(i==m)
                BB(m,2*m-1)=2;
                BB(m,2*m-2)=1;
        else
                BB(i,2*i-2)=1;
                BB(i,2*i-1)=2;
                BB(i,2*i)=1;
        end
    end
end
BB=BB./2;
