function [f,g]=sfun_XTM(W,DisR,MU_e,I0,L,m,nTau,NumElement)
global SigMa_XTM LogScale numThetan
global Tik penalty EM
%%===== Reconstruction discrete objective
%%===== L: intersection length matrix
%%===== M: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-M_i||^2, i=1..theta
mtol=prod(m);
lambda=1e-3;
L=reshape(full(L),numThetan,nTau+1,mtol);
%%%%% =================== Attenuation Matrix at beam energy
W=reshape(W,[1,1 mtol NumElement]);
MUe=squeeze(MU_e(:,1,1));
MU_XTM=squeeze(W)*MUe;
Mt=-log(DisR'./I0);%DisR';% 
Rdis=sum(bsxfun(@times,reshape(MU_XTM,[1,1,mtol]),L),3); %% Discrete case
%%%====================================
if(penalty)
    f=lambda*(norm(Tik*W))^2;
else
    f=0;
end
if(EM)
     Rdis = Rdis +1;
     Mt = Mt + 1;
     f_XTM=sum(sum((-log(Rdis).*Mt+Rdis),2),1);
     g_XTM =squeeze(sum(sum(bsxfun(@times,reshape(bsxfun(@times,1-Mt./Rdis,L),1,numThetan,nTau+1,mtol),MUe),2),3))';
else
     f_XTM=sum((Rdis(:)-Mt(:)).^2);
     g_XTM =2*squeeze(sum(sum(bsxfun(@times,reshape(bsxfun(@times,Rdis-Mt,L),1,numThetan,nTau+1,mtol),MUe),2),3))';
 end
g=g_XTM(:);
f=f+f_XTM;
if(penalty)
    Wdiff=W(1:end-1)-W(2:end);
    g=g+lambda.*(-2.*[0;Wdiff]+2.*[Wdiff;0]);
    % g=g+lambda*2*Tik'*(Tik*W);
end





