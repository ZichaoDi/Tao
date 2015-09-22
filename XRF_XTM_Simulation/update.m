function vh=update(vH,res_prob)
global N NumElement current_n
%== Full Weighting to interpolate
%== res_prob=0: upsample data; res_prob=1: upsample variable
j=find(N==current_n);
nH=sqrt(length(vH)/NumElement);
nh=nH*2-1;
% if(res_prob)
    vH=reshape(vH,nH,nH,NumElement);
    vH_b=zeros(nH+2,nH+2,NumElement);vH_b(2:nH+1,2:nH+1,:)=vH;
    vH=vH_b;
    vh=zeros(nh+2,nh+2,NumElement);
    wm=1/4;
    wl=1/8;
    wd=1/16;
    for i=2:nH+1
        for nj=2:nH+1
            vh(2*(i-1),2*(nj-1),:)=vH(i,nj,:)*wm;
            
            vh(2*(i-1)-1,2*(nj-1),:)=(vH(i-1,nj,:)+vH(i,nj,:))*wl;
            vh(2*(i-1)+1,2*(nj-1),:)=(vH(i,nj,:)+vH(i+1,nj,:))*wl;
            vh(2*(i-1),2*(nj-1)-1,:)=(vH(i,nj-1,:)+vH(i,nj,:))*wl;
            vh(2*(i-1),2*(nj-1)+1,:)=(vH(i,nj,:)+vH(i,nj+1,:))*wl;
            
            vh(2*(i-1)-1,2*(nj-1)-1,:)=(vH(i-1,nj,:)+vH(i,nj,:)+vH(i-1,nj-1,:)+vH(i,nj-1,:))*wd;
            vh(2*(i-1)+1,2*(nj-1)+1,:)=(vH(i+1,nj,:)+vH(i,nj,:)+vH(i+1,nj+1,:)+vH(i,nj+1,:))*wd;
            vh(2*(i-1)-1,2*(nj-1)+1,:)=(vH(i-1,nj,:)+vH(i,nj,:)+vH(i,nj+1,:)+vH(i-1,nj+1,:))*wd;
            vh(2*(i-1)+1,2*(nj-1)-1,:)=(vH(i+1,nj,:)+vH(i,nj,:)+vH(i,nj-1,:)+vH(i+1,nj-1,:))*wd;
        end
    end
    
% end
vh=vh(2:nh+1,2:nh+1,:);
% figure(11);clims=[0 max([vh(:);vH(:)])];subplot(1,2,1),imagesc(2*vh,clims);
% subplot(1,2,2),imagesc(vH,clims);
% pause;
vh=4*vh(:);
