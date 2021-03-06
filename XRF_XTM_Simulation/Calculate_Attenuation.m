function [ConstSub, InTens, OutTens, AttenuM, DW,f]=Calculate_Attenuation(W_rep,NumElement,L,GlobalInd,SelfInd,m,nTau,xrfData,M)
%%==== Given elemental map W and pre-calculate the beam and fluorescent attenuation coefficients
global NumSSDlet numChannel NoSelfAbsorption numThetan
global icycle
mtol=prod(m);
L=reshape(full(L),numThetan,nTau+1,mtol);
%%%%% ====================================================================
f=0;
InTens=ones(numThetan,nTau+1,mtol);
OutTens=ones(numThetan,nTau+1,mtol,NumElement);
for n=1:numThetan
    OutTens_d=ones(nTau+1,mtol,NumSSDlet,NumElement);
    sum_Tau=0;
    for i=1:nTau+1
        XRF_v=zeros(1,numChannel);
        if(ndims(W_rep)==4)
            W=reshape(W_rep(n,i,:,:),mtol,NumElement);
        else
            W=reshape(W_rep,mtol,NumElement);
        end
        counted_v=zeros(mtol,1);
        index=GlobalInd{n,i};
        if(~isempty(index))
            index_sub=sub2ind(m,index(:,2),index(:,1));
            msub=length(index_sub);
            for v_count=1:msub
                v=index_sub(v_count);
                scaleMU=10^(6);
                if(~isempty(SelfInd{n,i,v}{1}))
                    InTens(n,i,v)=exp(-sum(sum(W(SelfInd{n,i,v}{1},:).*SelfInd{n,i,v}{3}*scaleMU)));
                end
                if(~isempty(SelfInd{n,i,v}{2})& ~NoSelfAbsorption)
                    for d=1:NumSSDlet
                        if(~isempty(SelfInd{n,i,v}{2}{d}))
                        % temp_out=-sum(sum(bsxfun(@times,W(SelfInd{n,i,v}{2}{d},:),reshape(SelfInd{n,i,v}{4}{d},length(SelfInd{n,i,v}{2}{d}),NumElement,NumElement)),1),2);
                        OutTens_d(i,v,d,:)=exp(-sum(sum(bsxfun(@times,W(SelfInd{n,i,v}{2}{d},:),reshape(SelfInd{n,i,v}{4}{d}*scaleMU,length(SelfInd{n,i,v}{2}{d}),NumElement,NumElement)),1),2));
                        end
                    end
                OutTens(n,i,v,:)=sum(OutTens_d(i,v,:,:),3)/NumSSDlet;
                end
                XRF_v=XRF_v+L(n,i,v)*(InTens(n,i,v)*reshape(OutTens(n,i,v,:),1,NumElement).*W(v,:))*M;
            end
            sum_Tau=sum_Tau+(squeeze(xrfData(n,i,:))'-XRF_v)*(squeeze(xrfData(n,i,:))'-XRF_v)';
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
clear L_rep temp_v Mrep_P
