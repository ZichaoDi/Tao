%%==== Given elemental map W and pre-calculate the beam and fluorescent attenuation coefficients
global TempBeta Beta Joint frame
mtol=prod(m);
W=reshape(W,mtol,NumElement);
MU=zeros(prod(m),NumElement);
for i=1:NumElement
    temp=sum(reshape(W,prod(m),NumElement)*MU_e(:,:,i+1),2);
    temp=flipud(reshape(temp,m(1),m(2))');
    MU(:,i)=temp(:);
end
ConstSub=sparse(numThetan*(nTau+1)*numChannel,mtol*NumElement);
xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
%%%%% ====================================================================
    for n=1:numThetan
        theta=thetan(n)/180*pi;
        fprintf(1,'====== Angle Number  %d of %d: %d\n',n,numThetan,thetan(n));
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
            index=GlobalInd{ind_bt};
            if(~isempty(index))
                CurrentCellCenter=[(index(:,1)-1/2)*dz(1)-abs(omega(1)),(index(:,2)-1/2)*dz(2)-abs(omega(3))];
                currentInd=sub2ind(m,index(:,2),index(:,1));
                for j=1:length(currentInd)
                    InTens=exp(-sum(sum(W(sub2ind(m,index(j-1:-1:1,2),index(j-1:-1:1,1)),:).*kron(Lvec(j-1:-1:1),MU_e(:,1,1)'))));
                    in_after=setdiff(find(inpolygon(xc(:,1),xc(:,2),[CurrentCellCenter(j,1) SSDknot(1,1) SSDknot(NumSSDlet,1)],[CurrentCellCenter(j,2) SSDknot(1,2) SSDknot(NumSSDlet,2)])),currentInd); %% energy is not attenuated from the source point
                    if( ~NoSelfAbsorption)
                        OutTens=exp(-sum(MU(in_after,1:NumElement),1)./(length(in_after)+1)*prod(dz)*length(in_after));
                    end
                    ConstSub(sub2ind([numThetan*(nTau+1),numChannel],ind_bt*ones(1,numChannel),1:numChannel),sub2ind([mtol,NumElement],currentInd(j)*ones(1,NumElement),1:NumElement))=M'.*repmat(Lvec(j)*InTens*OutTens,numChannel,1);

                end
            end
        end
    end
%#############################################
XRF_v=ConstSub*W(:);
if(Joint==0)
    if(strcmp(frame,'EM'))
        thres=1;
        f=sum(-log(XRF_v+thres).*(xrfData(:)+thres)+XRF_v+thres);
    else
        f=sum((XRF_v-xrfData).^2);  
    end
elseif(Joint==1)
    MU_XTM=W*squeeze(MU_e(:,1,1));
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
end
