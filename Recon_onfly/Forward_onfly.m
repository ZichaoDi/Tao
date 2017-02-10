function [ConstSub, f]=Forward_onfly(x,xrfData,Mt, M)
%%==== Given elemental map x and pre-calculate the beam and fluorescent attenuation coefficients
global initialize TempBeta Beta Joint frame
global NoSelfAbsorption
global m numChannel NumElement nTau numThetan NumSSDlet
global SSDlet DetKnot0 SourceKnot0 thetan omega 
global MU_e dz mtol xbox ybox xc BeforeEmit 
global L in_after
%%%=========Start Simulation
x=reshape(x,mtol,NumElement);
MU=zeros(mtol,NumElement);
for i=1:NumElement
    temp=sum(x*MU_e(:,:,i+1),2);
    temp=flipud(reshape(temp,m(1),m(2))');
    MU(:,i)=temp(:);
end
ConstSub=sparse(numThetan*(nTau+1)*numChannel,mtol*NumElement);
%%%%% ====================================================================
fprintf(1,'====== Start Forward Mapping: %d angles %d beamlets\n',numThetan, nTau+1);
    for n=1:numThetan
        if(initialize)
            fprintf(1,'====== Angle Number  %d of %d: %d\n',n,numThetan,thetan(n));
        end
        theta=thetan(n)/180*pi;
        TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
        DetKnot=DetKnot0*TransMatrix;
        SourceKnot=SourceKnot0*TransMatrix;
        SSDknot=SSDlet*TransMatrix;
        %=================================================================
        CurrentCellCenter=[];
        CurrentInd=[];
        for i=1:nTau+1
            % Initialize
            BeforeEmit=1;
            %=================================================================
            [index,Lvec,linearInd]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta);
            ind_bt=(i-1)*numThetan+n;
            if(~isempty(index))
                CurrentCellCenter=[(index(:,1)-1/2)*dz(1)-abs(omega(1)),(index(:,2)-1/2)*dz(2)-abs(omega(3))];
                currentInd=sub2ind(m,index(:,2),index(:,1));
                if(initialize)
                    L(ind_bt,currentInd)=Lvec;
                end
                InTens=exp(-cumsum([0;sum(x(currentInd(1:end-1),:).*kron(Lvec(1:end-1),MU_e(:,1,1)'),2)]));
                OutTens=[];
                for j=1:length(currentInd)
                    if(initialize)
                        in_after{ind_bt,j}=int32(setdiff(find(inpolygon(xc(:,1),xc(:,2),[CurrentCellCenter(j,1) SSDknot(1,1) SSDknot(NumSSDlet,1)],[CurrentCellCenter(j,2) SSDknot(1,2) SSDknot(NumSSDlet,2)])),currentInd(j))); %% energy is not attenuated from the source point
                    end
                    if( ~NoSelfAbsorption)
                        OutTens(j,:)=exp(-sum(MU(in_after{ind_bt,j},1:NumElement),1)./(length(in_after{ind_bt,j})+1)*prod(dz)*length(in_after{ind_bt,j}));
                    end
                end
                ConstSub(ind_bt+([1:numChannel]-1)*(numThetan*(nTau+1)),reshape(repmat(currentInd,1,NumElement)'+(repmat([1:NumElement],length(currentInd),1)'-1).*mtol,length(currentInd)*NumElement,1))=...
                   reshape(repmat(M',[1,length(currentInd),1]),[numChannel,length(currentInd)*NumElement]).*repmat(reshape(repmat(Lvec.*InTens,1,NumElement)'.*OutTens',[length(currentInd)*NumElement,1])',numChannel,1);
            end
        end
    end
%#############################################
    XRF_v=ConstSub*x(:);
    MU_XTM=x*squeeze(MU_e(:,1,1));
    Rdis=reshape(L,numThetan*(nTau+1),mtol)*MU_XTM;
    if(strcmp(frame,'EM'))
        thres=1;
        f=sum(-log(XRF_v+thres).*(xrfData(:)+thres)+XRF_v+thres);
        f_XTM=sum(-log(Rdis+thres).*(Mt+thres)+Rdis+thres);
    elseif(strcmp(frame,'LS'))
        f=sum((XRF_v-xrfData(:)).^2);  
        f_XTM=sum((Rdis-Mt).^2);
    elseif(strcmp(frame,'mix'))
        thres=1;
        f=sum(-log(XRF_v+thres).*(xrfData(:)+thres)+XRF_v+thres);
        f_XTM=sum((Rdis-Mt).^2);
    end
    f = TempBeta*f + Beta*f_XTM;
