function [f,g]=sfun_AdiMat(W,xrfData,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau)

opts = admOptions('independents', [1]);
f=feval(@func_Tensor,W,xrfData,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau);
g=ones(1,length(W))/length(W)*admDiffRev(@func_Tensor,ones(size(W)),W,xrfData,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,opts);
g=g';
