function vh=update(vH,res_prob)
global N NumElement
%== Full Weighting to interpolate
%== res_prob=0: upsample data; res_prob=1: upsample variable
j=find(N==sqrt(size(vH,1)/NumElement));
% if(res_prob)
    vH=reshape(vH,N(j),N(j),NumElement);
    vH_b=zeros(N(j)+2,N(j)+2,NumElement);vH_b(2:N(j)+1,2:N(j)+1,:)=vH;
    vH=vH_b;
    vh=zeros(N(j-1)+2,N(j-1)+2,NumElement);
    for i=2:N(j)+1
        for nj=2:N(j)+1
            vh(2*(i-1),2*(nj-1),:)=vH(i,nj,:);
            
            vh(2*(i-1)-1,2*(nj-1),:)=(vH(i-1,nj,:)+vH(i,nj,:))/2;
            vh(2*(i-1)+1,2*(nj-1),:)=(vH(i,nj,:)+vH(i+1,nj,:))/2;
            vh(2*(i-1),2*(nj-1)-1,:)=(vH(i,nj-1,:)+vH(i,nj,:))/2;
            vh(2*(i-1),2*(nj-1)+1,:)=(vH(i,nj,:)+vH(i,nj+1,:))/2;
            
            vh(2*(i-1)-1,2*(nj-1)-1,:)=(vH(i-1,nj,:)+vH(i,nj,:)+vH(i-1,nj-1,:)+vH(i,nj-1,:))/4;
            vh(2*(i-1)+1,2*(nj-1)+1,:)=(vH(i+1,nj,:)+vH(i,nj,:)+vH(i+1,nj+1,:)+vH(i,nj+1,:))/4;
            vh(2*(i-1)-1,2*(nj-1)+1,:)=(vH(i-1,nj,:)+vH(i,nj,:)+vH(i,nj+1,:)+vH(i-1,nj+1,:))/4;
            vh(2*(i-1)+1,2*(nj-1)-1,:)=(vH(i+1,nj,:)+vH(i,nj,:)+vH(i,nj-1,:)+vH(i+1,nj-1,:))/4;
        end
    end
    
% end
vh=vh(2:N(j-1)+1,2:N(j-1)+1,:);
vh=2*vh(:);
