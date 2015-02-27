function vH=downdate(vh,res_prob)
global N current_n NumElement
%== Full Weighting to restrict
%== res_prob=0: downsample data; res_prob=1: downsample variable
j=find(N==current_n);
% if(res_prob)
    vh=reshape(vh,N(j),N(j),NumElement);
    vH=zeros(N(j+1),N(j+1),NumElement);
    vh_b=zeros(N(j)+2,N(j)+2,NumElement);vh_b(2:N(j)+1,2:N(j)+1,:)=vh;
    for i=1:N(j+1)
        for nj=1:N(j+1)
            vH(i,nj,:)=vh_b(2*i,2*nj,:)+...
                (vh_b(2*i,2*nj-1,:)+vh_b(2*i,2*nj+1,:)+vh_b(2*i+1,2*nj,:)+vh_b(2*i-1,2*nj,:))/2+...
                (vh_b(2*i-1,2*nj-1,:)+vh_b(2*i-1,2*nj+1,:)+vh_b(2*i+1,2*nj-1,:)+vh_b(2*i+1,2*nj+1,:))/4;           
        end
    end
    
% end
vH=vH(:)/2;
