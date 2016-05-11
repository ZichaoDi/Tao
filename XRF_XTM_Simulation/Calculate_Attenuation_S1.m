function [ConstSub, f]=Calculate_Attenuation(W_rep,NumElement,L,GlobalInd,SelfInd,m,nTau,xrfData,M)
%%==== Given elemental map W and pre-calculate the beam and fluorescent attenuation coefficients
global NumSSDlet numChannel NoSelfAbsorption numThetan
global icycle maxOut MU_e area_xrf EM I0
mtol=prod(m);
L=reshape(full(L),numThetan,nTau+1,mtol);
W=reshape(W_rep,mtol,NumElement);
MU=W*reshape(MU_e(:,1,:),NumElement,NumElement+1);
%%%%% ====================================================================
f=0;
totPro=numThetan*(nTau+1);
InTens=ones(numThetan,nTau+1,mtol);
OutTens=ones(numThetan,nTau+1,mtol,NumElement);
for n=1:numThetan
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
                    InTens(n,i,v)=exp(-sum(sum(W(SelfInd{n,i,v}{1},:).*SelfInd{n,i,v}{3})));
                end
                if( ~isempty(SelfInd{n,i,v}{2})& ~NoSelfAbsorption)
                    OutTens(n,i,v,:)=exp(-sum(MU(SelfInd{n,i,v}{2},2:NumElement+1),1)./(length(SelfInd{n,i,v}{2})+1)*area_xrf(n,i,v));
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
L_rep=reshape(L,numThetan,nTau+1,mtol);
temp_v=L_rep.*InTens;
Mrep_P=repmat(reshape(M,[1,1,1 NumElement numChannel]),[numThetan,nTau+1,mtol,1,1]);
ConstSub=bsxfun(@times,bsxfun(@times,temp_v,OutTens),Mrep_P); % part 3
clear L_rep temp_v Mrep_P OutTens InTens
ConstSub=reshape(permute(ConstSub,[1 2 5 3 4]),[numThetan*(nTau+1)*numChannel,mtol*NumElement]);
