function BB=restric_residule(n)
if(mod(n,2)==0)
    m=n/2;
else
    m=(n-1)/2;
end
BB=sparse(zeros(m,n));
for i=1:m
    if(i==m&mod(n,2)==0)
        BB(m,2*m-1)=1;
        BB(m,2*m)=2;
    else
        BB(i,2*i-1)=1;
        BB(i,2*i)=2;
        BB(i,2*i+1)=1;
   end
end