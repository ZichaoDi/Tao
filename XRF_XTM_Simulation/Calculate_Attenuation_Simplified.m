function [ConstSub, InTens, OutTens, AttenuM, DW,f]=Calculate_Attenuation(W_rep,NumElement,L,GlobalInd,SelfInd,m,nTau,xrfData,M)
%%==== Given elemental map W and pre-calculate the beam and fluorescent attenuation coefficients
global NumSSDlet numChannel NoSelfAbsorption numThetan
global icycle maxOut MU_e area_xrf EM scaleMU pix_inner
mtol=prod(m);
L=reshape(full(L),numThetan,nTau+1,mtol);
W=reshape(W_rep,mtol,NumElement);
MU=W*reshape(MU_e(:,1,:),NumElement,NumElement+1);
% if(icycle ==1)
%    MU(pix_inner,3)=max(MU(pix_inner,3));
%    MU(setdiff([1:mtol],pix_inner),3)=0;
% end
%%%%% ====================================================================
f=0;
InTens=ones(numThetan,nTau+1,mtol);
OutTens=ones(numThetan,nTau+1,mtol,NumElement);
for n=1:numThetan
    OutTens_d=ones(nTau+1,mtol,NumSSDlet,NumElement);
    sum_Tau=0;
    for i=1:nTau+1
        XRF_v=zeros(1,numChannel);
        counted_v=zeros(mtol,1);
        index=GlobalInd{n,i};
        if(~isempty(index))
            index_sub=sub2ind(m,index(:,2),index(:,1));
            msub=length(index_sub);
            for v_count=1:msub
                v=index_sub(v_count);
                if(~isempty(SelfInd{n,i,v}{1}))
                    InTens(n,i,v)=exp(-sum(sum(W(SelfInd{n,i,v}{1},:).*SelfInd{n,i,v}{3}*scaleMU)));
                end
                if( ~isempty(SelfInd{n,i,v}{2})& ~NoSelfAbsorption)
                    OutTens(n,i,v,:)=exp(-sum(MU(SelfInd{n,i,v}{2},2:NumElement+1),1)./(length(SelfInd{n,i,v}{2})+1)*area_xrf(n,i,v)*scaleMU);
                end
                XRF_v=XRF_v+L(n,i,v)*(InTens(n,i,v)*reshape(OutTens(n,i,v,:),1,NumElement).*W(v,:))*M;
            end
            if(EM)
                XRF_v=XRF_v+1;
                xrfData=xrfData+1;
                sum_Tau=sum_Tau+sum((-log(XRF_v).*reshape(xrfData(n,i,:),1,numChannel)+XRF_v),2);
            else
                sum_Tau=sum_Tau+(squeeze(xrfData(n,i,:))'-XRF_v)*(squeeze(xrfData(n,i,:))'-XRF_v)';     
            end
        end
    end
    f=f+sum_Tau;
end
AttenuM=bsxfun(@times,InTens,OutTens);
DW=bsxfun(@times,AttenuM,reshape(W,[1,1,mtol,NumElement]));
L_rep=reshape(L,numThetan,nTau+1,mtol);
temp_v=L_rep.*InTens;
Mrep_P=repmat(reshape(M,[1,1,1 NumElement numChannel]),[numThetan,nTau+1,mtol,1,1]);
ConstSub=bsxfun(@times,bsxfun(@times,temp_v,OutTens),Mrep_P); % part 3
% O2=squeeze(sum(sum(OutTens,1),2));O22=map1D(O2,[0.01 1]);
% pert = 1e-1;
% scale = 3;
% xnew   = map1D((W+pert)./(map1D(O2,[0.1 1]).^scale),[min(W(:)),max(W(:))]);
% nrow=3;
% figure, for ii=1:3, subplot(nrow,3,ii);imagesc(reshape(W(:,ii),m(1),m(1)));colorbar;
% subplot(nrow,3,3+ii),imagesc(reshape(O2(:,ii),m(1),m(1)));colorbar;
% subplot(nrow,3,6+ii);imagesc(reshape(xnew(:,ii),m(1),m(1)));colorbar;
% if(nrow==4)
%     subplot(nrow,3,9+ii),imagesc(reshape(MU(:,ii+1),m(1),m(1)));colorbar;
% end
% end
% drawnow;
clear L_rep temp_v Mrep_P
