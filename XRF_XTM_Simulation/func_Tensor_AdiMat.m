function f=func_Tensor_AdiMat(W,xrfData,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau)
global NumSSDlet numChannel NoSelfAbsorption mtol
f=0;
W=reshape(W,mtol,NumElement);
L=reshape(L,length(thetan),nTau+1,mtol);
%%%%% ====================================================================
for n=1:length(thetan)
    sum_Tau=0;
    for i=1:nTau+1
        XRF_v=zeros(1,numChannel);
        index=GlobalInd{n,i};
        if(~isempty(index))
            index_sub=sub2ind(m,index(:,2),index(:,1));
            msub=length(index_sub);
            for v_count=1:msub
                v=index_sub(v_count);
                if(isempty(SelfInd{n,i,v}{1}))
                    InTens=1;
                else
                    InTens=exp(-sum(sum(W(SelfInd{n,i,v}{1},:).*SelfInd{n,i,v}{3})));
                end
                
                if(isempty(SelfInd{n,i,v}{2})|NoSelfAbsorption)
                    OutTens_d=ones(1,NumElement);
                else
                    OutTens_d=zeros(1,1,NumElement);
                    for d=1:NumSSDlet
                        if(~isempty(SelfInd{n,i,v}{2}{d}))
                            Tmmp=repmat(W(SelfInd{n,i,v}{2}{d},:),[1,1,NumElement]).*reshape(SelfInd{n,i,v}{4}{d},length(SelfInd{n,i,v}{2}{d}),NumElement,NumElement);
                        else
                            Tmmp=0;
                        end
                        OutTens_d=OutTens_d+exp(-sum(sum(Tmmp,1),2));
                        %                  OutTens_d=OutTens_d+exp(-sum(sum(bsxfun(@times,W(SelfInd{n,i,v}{2}{d},:),reshape(SelfInd{n,i,v}{4}{d},length(SelfInd{n,i,v}{2}{d}),NumElement,NumElement)),1),2));
                    end
                    OutTens_d=OutTens_d/NumSSDlet;
                end
                XRF_v=XRF_v+L(n,i,v)*(InTens*reshape(OutTens_d,1,NumElement).*W(v,:))*M;
            end
%             sum_Tau
            sum_Tau=sum_Tau+(xrfData{n,i}-XRF_v)*(xrfData{n,i}-XRF_v)';
        end
    end
    f=f+sum_Tau;
end
