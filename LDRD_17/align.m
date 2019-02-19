function S = align(x,W)
N=size(W,1);
clims=[min(W(:)),max(W(:))];
x=map1D(reshape(x,N,N),clims);
tf=imregcorr(x,W,'translation');
Rfixed = imref2d(size(W));
S=imwarp(x,tf,'OutputView',Rfixed);


