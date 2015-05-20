function vH=downdate(vh,res_prob)
global N current_n NumElement numThetan LevelScale
%== Full Weighting to restrict
%== res_prob => 0: downsample variable; 1: downsample gradient; 3: downsample data
j=find(N==current_n);
if(res_prob==3) %% vh: numThetan X nTau+1 X numChannel
    
    for n=1:numThetan
        if(ndims(vh)==2)
            vH(:,n)=restric_residule(size(vh,1))*vh(:,n);
        elseif(ndims(vh)==3)
            vH(n,:,:)=restric_residule(size(vh,2))*vh(n,:,:);
        end
    end
else
     NumElement=1;
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
    if(res_prob==0)
    vH=vH(:)/2;
    elseif(res_prob==1)
        vH=vH(:)/2;
    end
end


