function [f,g]=sfun_AdiMat(W,xrfData,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau)

opts = admOptions('independents', [1]);
f=feval(@func_Tensor_AdiMat,W,xrfData,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau);
if(NumElement==1)
g=ones(1,length(W))/length(W)*admDiffRev(@func_Tensor_AdiMat,ones(size(W)),W,xrfData,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,opts);
% g=admDiffFor(@func_Tensor_AdiMat,1,W,xrfData,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,opts);
else
g=admDiffFD(@func_Tensor_AdiMat,1,W,xrfData,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,opts);
end
g=g';

