function B=downdate_L(N,level)
%== Full Weighting to restrict
  %  vh=reshape(vh,N(level-1),N(level-1));
  %  vH=zeros(N(level),N(level));
  %  vh_b=zeros(N(level-1)+2,N(level-1)+2);vh_b(2:N(level-1)+1,2:N(level-1)+1)=vh;
  %  for i=1:N(level)
  %      for nj=1:N(level)
  %          vH(i,nj)=vh_b(2*i,2*nj)+...
  %              (vh_b(2*i,2*nj-1)+vh_b(2*i,2*nj+1)+vh_b(2*i+1,2*nj)+vh_b(2*i-1,2*nj))/2+...
  %              (vh_b(2*i-1,2*nj-1)+vh_b(2*i-1,2*nj+1)+vh_b(2*i+1,2*nj-1)+vh_b(2*i+1,2*nj+1))/4;
  %      end
  %  end
  %  vH=vH(:)/4;
w=[1/4, 1/8, 1/16]*1;
mh=[N(level-1)+2,N(level-1)+2];
mH=[N(level),N(level)];
B=zeros(prod(mH),prod(mh));
hInd=[1:N(level-1)+2]';
sx0=0*hInd+1; sx1=hInd; 
boarder=[sx0,sx1;sx1,sx0;sx0+hInd(end)-1,sx1;sx1,sx0+hInd(end)-1];
inner=setdiff(1:prod(mh),sub2ind(mh,boarder(:,1),boarder(:,2)));
for i=1: N(level)
    for j=1:N(level)
        B(sub2ind(mH,i,j),sub2ind(mh,2*i,2*j))=w(1);
        %% = vertical & horizontal
        B(sub2ind(mH,i,j),sub2ind(mh,2*i,2*j-1))=w(2);
        B(sub2ind(mH,i,j),sub2ind(mh,2*i,2*j+1))=w(2);
        B(sub2ind(mH,i,j),sub2ind(mh,2*i+1,2*j))=w(2);
        B(sub2ind(mH,i,j),sub2ind(mh,2*i-1,2*j))=w(2);
        %% = diagonal
        B(sub2ind(mH,i,j),sub2ind(mh,2*i-1,2*j-1))=w(3);
        B(sub2ind(mH,i,j),sub2ind(mh,2*i-1,2*j+1))=w(3);
        B(sub2ind(mH,i,j),sub2ind(mh,2*i+1,2*j-1))=w(3);
        B(sub2ind(mH,i,j),sub2ind(mh,2*i+1,2*j+1))=w(3);
    end
end
B=B(:,inner);
        


