function [f,g]=sfun_XRF(W,xrfData,MU,M,NumElement,numChannel,Ltol,GlobalInd,thetan,m,nTau)
f=0;
W=reshape(W,m(1),m(2),NumElement);
g=zeros(m(1),m(2),NumElement);
for n=1:length(thetan)
    sum_Tau=0;
    for i=1:nTau+1
        index=GlobalInd{n,i};
        if(~isempty(index))
            L=Ltol{n,i};
            RM=cell(m(1),m(2));
            xrfSub=zeros(1,numChannel);
            TempSub=zeros(NumElement,size(index,1),numChannel);
            for j=1:size(index,1)
                
                if(j==1)
                    I_incident=1;
                    temp_sum=0;
                else
                    
                    temp_sum=temp_sum+L(index(j-1,2),index(j-1,1))*MU(index(j-1,2),index(j-1,1));
                    I_incident=exp(-temp_sum);
                end
                Wsub=reshape(W(index(j,2),index(j,1),:),[NumElement,1]);
                RM{index(j,2),index(j,1)}=L(index(j,2),index(j,1))*I_incident*(Wsub'*M);
                TempSub(:,j,:)=(L(index(j,2),index(j,1))*I_incident).*M;
                     xrfSub=xrfSub+RM{index(j,2),index(j,1)};
            end
            for j=1:size(index,1)
                temp=zeros(size(g(index(j,2),index(j,1),:)));
                temp(1,1,:)=2*reshape(TempSub(:,j,:),NumElement,numChannel)*(xrfSub-xrfData{n,i})';
                g(index(j,2),index(j,1),:)=g(index(j,2),index(j,1),:)+temp;
               
            end
            sum_Tau=sum_Tau+(xrfData{n,i}-xrfSub)*(xrfData{n,i}-xrfSub)';
        end
    end
    f=f+sum_Tau;
end
g=g(:);
