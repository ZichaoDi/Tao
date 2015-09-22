%----------------------------------------------------------------------
% Perform "global" computations for optimization at a
% specified discretization (corresponding to n)

global N xrf_level xtm_level L_level GI_level SI_level SigmaR SigmaT m_level nTau_level W_level
global SigMa_XTM SigMa_XRF L GlobalInd SelfInd m nTau XRF DisR
global W0 thetan Joint MUe WS ReconAttenu

level=find(N==current_n);
GlobalInd=GI_level{level};
SelfInd=SI_level{level};
m=m_level(level,:);
nTau=nTau_level(level);
if(Joint==-1)
    L=L_level{level};
else
    L=reshape(full(L_level{level}),length(thetan),nTau+1,prod(m));
end
SigMa_XTM=SigmaT{level};
SigMa_XRF=SigmaR{level};
if(level>1)
    W0=W_level{level};
    W0=W0(:);
else
    W0=WS(:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(level==1 | init_level==1)
    XRF=xrf_level{1};
    DisR=xtm_level{1};
    %%%----------------------------------------------------------------------
elseif(level~=1 & init_level == 0)
    %     if(level~=1 & init_level == 0)
        
    [Fv2 ,Gv2 , shift_yH, shift_yTH] = sfun(current_v);
    if(level==2)
        LevelScale(level-1)=1;%abs(Fv2/Fv);
    else
        LevelScale(level-1)=1;%abs(Fv2/Fv);
    end
    %     end
    if(Joint==1)
        XRF=XRF-shift_y+shift_yH;
        xrf_level{level}=XRF;
        DisR=DisR-shift_yT'+shift_yTH';
        xtm_level{level}=DisR;
    elseif(Joint==0)
        XRF=XRF-shift_y+shift_yH;
        xrf_level{level}=XRF;
    elseif(Joint==-1)
         DisR=DisR-shift_yT'+shift_yTH';
        xtm_level{level}=DisR;
    end
end

