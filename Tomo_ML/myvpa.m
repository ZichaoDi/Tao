function newA=myvpa(A)
newA=zeros(size(A));
for i=1:length(A)
if(A(i)<1e-10)
    newA(i)=0;
else
    newA(i)=A(i);
end
end