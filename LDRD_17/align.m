function S = align(x,W)
[n,m]=size(W);
clims=[min(W(:)),max(W(:))];
x=map1D(reshape(x,n,m),clims);
tf=imregcorr(x,W,'translation');
Rfixed = imref2d(size(W));
S=imwarp(x,tf,'OutputView',Rfixed);


