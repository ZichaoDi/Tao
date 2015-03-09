%----- Define Geometry for a given grid, thetan and experimental equipment
Define_Detector_Beam_Gaussian; %% provide the beam source and Detectorlet
thetan=linspace(0,180,30);%mod(thetan+360,360);% Projection Angles, has to be positive.
thetan=thetan(1:numThetan);
%%---------------------------------------------------------------------------
ID=cell(numThetan,nTau+1);
LD=cell(numThetan,nTau+1);
LA=cell(numThetan,nTau+1,m(1),m(2),NumSSDlet);
LI=cell(numThetan,nTau+1,m(1),m(2),NumSSDlet);
for n=1:numThetan
    theta=thetan(n)/180*pi;
fprintf(1,'====== Angle Number  %d of %d: %d\n',n,numThetan,thetan(n));
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    DetKnot=DetKnot0*TransMatrix;
    SourceKnot=SourceKnot0*TransMatrix;
    SSDknot=SSDlet*TransMatrix;
    %%---------------------------------------------------------------------------
    for i=1:nTau+1 %%%%%%%%%================================================================
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        BeforeEmit=1;
        %====================================
%         figure('name',sprintf('XRF at beam %d with angle %d',i,thetan(n)));
%             subplot(1,2,1);
%             xc = getNodalGrid(omega,[m(2) m(1)]);
%             plotGrid(xc,omega,[m(2) m(1)]);
%             hold on;
%             plot(DetKnot(:,1),DetKnot(:,2),'k+-',SourceKnot(:,1),SourceKnot(:,2),'m+-',SSDknot(:,1),SSDknot(:,2),'g+-','LineWidth',0.5)
%             axis equal
%             set(gcf,'Units','normalized')
%             set(gca,'Units','normalized')
%             ax = axis;
%             ap = get(gca,'Position');
%             xp = ([SourceKnot(i,1),DetKnot(i,1)]-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
%             yp = ([SourceKnot(i,2),DetKnot(i,2)]-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
%             ah=annotation('arrow',xp,yp,'Color','r','LineStyle','--');
        
        %=================================================================
        [index,Lvec]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta);
        ID{n,i}=index;
        LD{n,i}=Lvec;
        %%%%%%%%================================================================
        for j=1:size(index,1)
            CurrentCellCenter=[(index(j,1)-1/2)*dz(1)-abs(omega(1)),(index(j,2)-1/2)*dz(2)-abs(omega(3))];
            BeforeEmit=0;
            for SSDi=1:NumSSDlet
                beta=angle(SSDknot(SSDi,1)-CurrentCellCenter(1)+(SSDknot(SSDi,2)-CurrentCellCenter(2))*sqrt(-1));
                if(beta>=0 & beta<=pi/2)
                    xbox=[CurrentCellCenter(1) CurrentCellCenter(1) omega(2) omega(2) CurrentCellCenter(1)];
                    ybox=[CurrentCellCenter(2) omega(4) omega(4) CurrentCellCenter(2) CurrentCellCenter(2)];
                elseif(beta >pi/2 & beta<=pi)
                    xbox=[omega(1) omega(1) CurrentCellCenter(1) CurrentCellCenter(1) omega(1)];
                    ybox=[ CurrentCellCenter(2) omega(4) omega(4) CurrentCellCenter(2) CurrentCellCenter(2)];
                elseif(beta >=-pi/2 & beta<0)
                    xbox=[CurrentCellCenter(1) CurrentCellCenter(1) omega(2) omega(2) CurrentCellCenter(1)];
                    ybox=[omega(3) CurrentCellCenter(2) CurrentCellCenter(2) omega(3) omega(3)];
                elseif(beta>=-pi & beta<-pi/2)
                    xbox=[omega(1) omega(1) CurrentCellCenter(1) CurrentCellCenter(1) omega(1)];
                    ybox=[omega(3) CurrentCellCenter(2) CurrentCellCenter(2) omega(3) omega(3)];
                end
                if(beta<0)
                    beta=beta+2*pi;
                end
                [index_after,Lvec_after,linearInd_after]=IntersectionSet(CurrentCellCenter,SSDknot(SSDi,:),xbox,ybox,beta);
                [index_after,otherInd]=setdiff(index_after,index(j,:),'rows');
                Lvec_after=Lvec_after(otherInd');
                LinearInd=sub2ind(m,index_after(:,2),index_after(:,1));
                LA{n,i,index(j,2),index(j,1),SSDi}=Lvec_after;
                LI{n,i,index(j,2),index(j,1),SSDi}=LinearInd;
                
                
            end %% End loop for existing fluorescence energy from current pixel
            %% ====================================================================
        end
    end
end




