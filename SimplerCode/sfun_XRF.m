function [f,g]=sfun_XRF(W,xrfData,I,M,Energy,EnergyChannel,Ltol,GlobalInd,thetan,m,nTau)
f=0;
NumElement=size(M,1);
W=reshape(W,m(1),m(2),NumElement);
g=zeros(m(1),m(2),NumElement);
for n=1:length(thetan)
    sum_Tau=0;
    for i=1:nTau+1
        index=GlobalInd{n,i};
        if(~isempty(index))
            L=Ltol{n,i};
            RM=cell(m(1),m(2));
            xrfSub=zeros(size(EnergyChannel));
            TempSub=zeros(NumElement,size(index,1));
            for j=1:size(index,1)
                
                if(j==1)
                    I_incident=1;
                    temp_sum=0;
                else
                    
                    temp_sum=temp_sum+L(index(j-1,2),index(j-1,1))*I(index(j-1,2),index(j-1,1));
                    I_incident=exp(-temp_sum);
                end
                Wsub=reshape(W(index(j,2),index(j,1),:),size(M));
                RM{index(j,2),index(j,1)}=L(index(j,2),index(j,1))*I_incident*(Wsub.*M);
                TempSub(:,j)=(L(index(j,2),index(j,1))*I_incident).*M;
                ExistInd=find(Wsub~=0);
                for tsub=1:length(ExistInd)
                    DetInd=find(Energy(ExistInd(tsub))==EnergyChannel);
                    
                    xrfSub(DetInd)=xrfSub(DetInd)+RM{index(j,2),index(j,1)}(ExistInd(tsub));
                end
                
                
                
            end
            p1=GaussianFit1(EnergyChannel,xrfSub);
            for j=1:size(index,1)
                temp=zeros(size(g(index(j,2),index(j,1),:)));
                temp(1,1,:)=2.*(p1-xrfData{n,i}).*TempSub(:,j);
                g(index(j,2),index(j,1),:)=g(index(j,2),index(j,1),:)+temp;
               
            end
         
            sum_Tau=sum_Tau+(xrfData{n,i}-p1)'*(xrfData{n,i}-p1);
        end
    end
    f=f+sum_Tau;
end
g=g(:);

