function [index,Lvec,linearInd]=IntersectionSet(Source,Detector,xbox,ybox,theta)
global x y omega m plotTravel BeforeEmit dz  
global  fig2 fig5 finalfig
[Ax, Ay] = polyxpoly([Source(1),Detector(1)],[Source(2),Detector(2)], xbox, ybox);
if(isempty(Ax) | (length(Ax)==1 & length(Ay)==1 ))
    %fprintf('no intersection \n')
    index=[];
    Lvec=[];
    linearInd=[];
else
     A=unique([Ax,Ay],'rows');Ax=A(:,1);Ay=A(:,2);
    if(theta==pi/2)
        Q=[repmat(Ax(1),size(y')),y'];
    elseif(theta==0 | theta==2*pi)
        Q=[x',repmat(Ay(1),size(x'))];
    elseif(theta==pi)
        Q=[x(end:-1:1)',repmat(Ay(1),size(x'))];
        
    elseif(theta==3*pi/2)
        Q=[repmat(Ax(1),size(y')),y(end:-1:1)'];
        
    else
        Q=[[x', (Ay(2)-Ay(1))/(Ax(2)-Ax(1)).*x'+(Ay(1)*Ax(2)-Ay(2)*Ax(1))/(Ax(2)-Ax(1))];...
            [(y'-(Ay(1)*Ax(2)-Ay(2)*Ax(1))/(Ax(2)-Ax(1)))./((Ay(2)-Ay(1))/(Ax(2)-Ax(1))),y']];
    end
    
    indx=find(Q(:,1)-xbox(1)<-1e-6 |Q(:,1)-xbox(3)>1e-6); indy=find(Q(:,2)-ybox(1)<-1e-6 |Q(:,2)-ybox(2)>1e-6);
    Q=setdiff(Q,Q([indx;indy],:),'rows');
    Q=unique(Q,'rows');
    if(BeforeEmit)
        dis=sqrt(sum(bsxfun(@minus,Q,Source).^2,2)); 
        [~,InterOrder]=sort(dis);
        Q=Q(InterOrder,:);
    else
        Q=unique(floor([Q;A].*1e4)./1e4,'rows');
    end
    Lvec=sqrt(bsxfun(@minus,Q(2:end,1),Q(1:end-1,1)).^2+bsxfun(@minus,Q(2:end,2),Q(1:end-1,2)).^2);
    %%%%%%%%%================================================================
    QC=(Q(1:end-1,:)+Q(2:end,:))/2;
    index=floor([(QC(:,1)-omega(1))/dz(1)+1, (QC(:,2)-omega(3))/dz(2)+1]);
    if(~BeforeEmit)
    index=index(index(:,1)>0 & index(:,1)<=m(2)& index(:,2)<=m(1) & index(:,2)>0,:);
    [~,subInd]=unique(index,'rows');
    index=index(sort(subInd),:);
    end
    linearInd=sub2ind(m,index(:,2),index(:,1));
    %%%%%%%%%================================================================
    if plotTravel
        if(BeforeEmit)
            
            figure(finalfig)
            subplot(1,2,1);
            set(fig2,'visible','off');
            drawnow;
            if(~isempty(index))
                fig2=plot((index(:,1)-1/2)*dz-abs(omega(1)),(index(:,2)-1/2)*dz-abs(omega(3)),'bo',Q(:,1),Q(:,2),'g-');
%         pause;        
            end
        else
            figure(finalfig)
            subplot(1,2,1);
            set(fig5,'visible','off');
            drawnow;
            if(~isempty(index))
                fig5=plot((index(:,1)-1/2)*dz-abs(omega(1)),(index(:,2)-1/2)*dz-abs(omega(3)),'bo',Q(:,1),Q(:,2),'g-'); 
            end
            pause;
        end
    end
    %%%%%%%%%================================================================
    
end