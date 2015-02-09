function sfun=global_setup (n)
global  xrf_level xtm_level L_level GI_level SI_level SigmaR SigmaT m_level nTau_level
global SigMa_XTM SigMa_XRF N
global NumElement MU_e I0 M thetan W_level W0
%----------------------------------------------------------------------
% Perform "global" computations for optimization at a 
% specified discretization (corresponding to n)
level=find(N==n);
XRF=xrf_level{level};
DisR=xtm_level{level};
L=L_level{level};
GlobalInd=GI_level{level};
SelfInd=SI_level{level};
m=m_level(level,:);
nTau=nTau_level(level);
SigMa_XTM=SigmaT{level};
SigMa_XRF=SigmaR{level};
W0=W_level{level};
W0=W0(:);
%----------------------------------------------------------------------
sfun=@(W)sfun_Tensor_Joint(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0);
