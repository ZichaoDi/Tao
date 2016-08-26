
global low up penalty
global W0 current_n
global fctn_f err0 fiter nit maxiter

%%%----------------------------Initialize dependent variables
level=1;
current_n=N(level);
W= W_level{level};
DisR=xtm_level{level};
L=L_level{level};
GlobalInd=GI_level{level};
m=m_level(level,:);
nTau=nTau_level(level);
%%%----------------------------------------------------------------------
fctn=@(W)sfun_radon(W,xtm_level{level},I0,L_level{level});
low=0*ones(current_n^2,1);
up=inf*ones(current_n^2,1);
maxiter=100;
x0=zeros(current_n^2,1);
W0=W_level{level};
err0=norm(W0-x0);
if(bounds)
[xstar,f,g,ierror] = tnbc (x0,fctn,low,up);
else
[xstar,f,g,ierror] = tn (x0,fctn);
end
doplot(1,xstar,W0,x0);
