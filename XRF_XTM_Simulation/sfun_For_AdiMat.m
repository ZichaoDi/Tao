function [f,g]=sfun_For_AdiMat(W,xrfData,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,LocalInd,L_after,thetan,m,nTau)

opts = admOptions('independents', [1]);
% tic;
f=feval(@func_for,W,xrfData,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,LocalInd,L_after,thetan,m,nTau);
% toc; 
% return;

% g=ones(1,length(W))/length(W)*admDiffRev(@func_Tensor_AdiMat,ones(size(W)),W,xrfData,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,opts);
% g=admDiffFor(@func_Tensor_AdiMat,1,W,xrfData,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,opts);
g=admDiffFD(@func_for,1,W,xrfData,MU_e,M,NumElement,numChannel,Ltol,GlobalInd,LocalInd,L_after,thetan,m,nTau,opts);
g=g';

