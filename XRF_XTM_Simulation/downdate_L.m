function vH=downdate(vh,level)
global N
%== Full Weighting to restrict
    vh=reshape(vh,N(level-1),N(level-1));
    vH=zeros(N(level),N(level));
    vh_b=zeros(N(level-1)+2,N(level-1)+2);vh_b(2:N(level-1)+1,2:N(level-1)+1)=vh;
    for i=1:N(level)
        for nj=1:N(level)
            vH(i,nj)=vh_b(2*i,2*nj)+...
                (vh_b(2*i,2*nj-1)+vh_b(2*i,2*nj+1)+vh_b(2*i+1,2*nj)+vh_b(2*i-1,2*nj))/2+...
                (vh_b(2*i-1,2*nj-1)+vh_b(2*i-1,2*nj+1)+vh_b(2*i+1,2*nj-1)+vh_b(2*i+1,2*nj+1))/4;
        end
    end
    vH=vH(:)/4;


