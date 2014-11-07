function f=func_Tensor(W,xrfData,M,NumElement,L,SelfInd,thetan,nTau)
global NumSSDlet numChannel NoSelfAbsorption
global mtol
f=0;
W=reshape(W,mtol,NumElement);
L=reshape(L,length(thetan),nTau+1,mtol);
%%%%% ====================================================================
for n=1:length(thetan)
    sum_Tau=0;
    for i=1:nTau+1
        XRF_v=zeros(1,numChannel);
         for v=1:mtol
            if(isempty(SelfInd{n,i,v}{1}))
                InTens=1;
            else                
                InTens=exp(-sum(sum(W(SelfInd{n,i,v}{1},:).*SelfInd{n,i,v}{3})));
            end
            
                if(isempty(SelfInd{n,i,v}{2})|NoSelfAbsorption)
                    OutTens_d=ones(1,NumElement);
                else
                OutTens_d=0;
                for d=1:NumSSDlet
                    OutTens_d=OutTens_d+exp(-sum(sum(bsxfun(@times,W(SelfInd{n,i,v}{2}{d},:),reshape(SelfInd{n,i,v}{4}{d},length(SelfInd{n,i,v}{2}{d}),NumElement,NumElement)),1),2));
                end
                OutTens_d=OutTens_d/NumSSDlet;
                end
                XRF_v=XRF_v+L(n,i,v)*(InTens*reshape(OutTens_d,1,NumElement).*W(v,:))*M;
        end
        sum_Tau=sum_Tau+(xrfData{n,i}-XRF_v)*(xrfData{n,i}-XRF_v)';
    end
    f=f+sum_Tau;
end
