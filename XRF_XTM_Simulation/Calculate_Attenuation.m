function [InTens, OutTens, AttenuM, DW]=Calculate_Attenuation(W_rep,NumElement,L,GlobalInd,SelfInd,m,nTau)
%%==== Given elemental map W and pre-calculate the beam and fluorescent attenuation coefficients
global NumSSDlet numChannel NoSelfAbsorption numThetan
mtol=prod(m);
L=reshape(L,numThetan,nTau+1,mtol);
%%%%% ====================================================================
InTens=ones(numThetan,nTau+1,mtol);
OutTens=ones(numThetan,nTau+1,mtol,NumElement);
for n=1:numThetan
    OutTens_d=ones(nTau+1,mtol,NumSSDlet,NumElement);
    for i=1:nTau+1
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
                if(~isempty(SelfInd{n,i,v}{1}))
                    InTens(n,i,v)=exp(-sum(sum(W(SelfInd{n,i,v}{1},:).*SelfInd{n,i,v}{3})));
                end
                if(~isempty(SelfInd{n,i,v}{2})& ~NoSelfAbsorption)
                    for d=1:NumSSDlet
                        if(~isempty(SelfInd{n,i,v}{2}{d}))
                        OutTens_d(i,v,d,:)=exp(-sum(sum(bsxfun(@times,W(SelfInd{n,i,v}{2}{d},:),reshape(SelfInd{n,i,v}{4}{d},length(SelfInd{n,i,v}{2}{d}),NumElement,NumElement)),1),2));
                        end
                    end
                OutTens(n,i,v,:)=sum(OutTens_d(i,v,:,:),3)/NumSSDlet;
                end
            end
        end
    end
end
AttenuM=bsxfun(@times,InTens,OutTens);
DW=bsxfun(@times,AttenuM,reshape(W,[1,1,mtol,NumElement]));
