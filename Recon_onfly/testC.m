function ConstSub=testC(N1,N2,N3,N4);
rng('default');
L=rand(N1,N2);
I=rand(N1,N2);
O=rand(N1,N2,N3);
M=rand(N3,N4); 

ConstSub=bsxfun(@times,reshape(bsxfun(@times,full(L.*I),O),[N1,1,N2,N3]),reshape(M',[1,N4,1, N3])); 

