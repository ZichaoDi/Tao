%----------------------------------------------------------------------
% Perform "global" computations for optimization at a
% specified discretization (corresponding to n)

global N xrf_level xtm_level L_level GI_level SI_level SigmaR SigmaT m_level nTau_level W_level
global SigMa_XTM SigMa_XRF L GlobalInd SelfInd m nTau XRF DisR
global W0 thetan Joint

level=find(N==current_n);
GlobalInd=GI_level{level};
SelfInd=SI_level{level};
m=m_level(level,:);
nTau=nTau_level(level);
L=reshape(L_level{level},length(thetan),nTau+1,prod(m));
SigMa_XTM=SigmaT{level};
SigMa_XRF=SigmaR{level};
W0=W_level{level};
W0=W0(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(level==1 | init_level==1)
    XRF=xrf_level{level};
    DisR=xtm_level{level};
%%%%%----------------------------------------------------------------------
elseif(level~=1 & init_level == 0)
% if(level~=1 & init_level == 0)
    [~ ,Gv2 ,shift_yH, shift_yTH] = sfun(current_v);
% end
%     XRF=XRF-shift_y+shift_yH;
    xrf_level{level}=XRF;
    if(Joint==1)
%     DisR=DisR-shift_yT'+shift_yTH';
    xtm_level{level}=DisR;
    end
end

