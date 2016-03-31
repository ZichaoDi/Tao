%%% Simulate XRF of a given fixed object with rotating predifined detector and beam
%%% Travelling of fluorescence photon is approximated as the area covered by solid angle 
global plotSpecSingle BeforeEmit plotTravel SSDlet area_xrf
global fig2  fig5 finalfig EmptyBeam RealBeam
global LogScale Tol
plotSpecSingle=0;
more off;
Define_Detector_Beam_Gaussian; %% provide the beam source and Detectorlet
DefineObject_Gaussian; % Produce W, MU_XTM
%%%%%%%==============================================================
if plotTravel
    fig2=[];  fig5=[];
end
%%%=========Start Simulation
mtol=prod(m);
eX=ones(m(1),1);
eY=ones(m(2),1);
if(synthetic)
    XRF=zeros(numThetan,nTau+1,numChannel);
else
    XRF_raw=zeros(numThetan,nTau+1,numChannel_raw);
    XRF_decom=zeros(numThetan,nTau+1,numChannel_decom);
end
DisR=zeros(nTau+1,numThetan);
L=zeros(numThetan,nTau+1,m(1),m(2));%zeros(numThetan*(nTau+1)*prod(m),1);
GlobalInd=cell(numThetan,nTau+1);
area_xrf=zeros(numThetan,nTau+1,mtol);
AS=cell(6,1);
SelfInd=repmat({AS},[numThetan,nTau+1,mtol]);
clear AS
EmptyBeam=[];
RealBeam=[];
if(synthetic)
    fprintf(1,'====== Fluorescence Detector Resolution is %d\n',numChannel);
else
    fprintf(1,'====== Fluorescence Detector Resolution is %d\n',numChannel_raw);
end
for n=1:numThetan
    theta=thetan(n)/180*pi;
    fprintf(1,'====== Angle Number  %d of %d: %d\n',n,numThetan,thetan(n));
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    DetKnot=DetKnot0*TransMatrix;
    SourceKnot=SourceKnot0*TransMatrix;
    SSDknot=SSDlet*TransMatrix;
    %{ 
    % ============================================================
             fig=figure('name',['grids & beam line at Angle',num2str(theta)]);
             % plotGrid(xc,omega,[m(2) m(1)]); hold on;
             for i=1:nTau+1
                 fig3=plot([SourceKnot(i,1),DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)],'r.-');hold on;
             end
             imagesc((x(1:end-1)+x(2:end))/2,(y(1:end-1)+y(2:end))/2,sum(W,3))
             pause(1);
    % ===================================================
    %} 
    
    Rdis=I0*ones(nTau+1,1);
    for i=1:nTau+1 %%%%%%%%%================================================================
        % Initialize
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        if(synthetic)
            XRF(n,i,:) = zeros(numChannel,1);
        else
            XRF_raw(n,i,:) = zeros(numChannel_raw,1);
            XRF_decom(n,i,:) = zeros(numChannel_decom,1);
        end
        BeforeEmit=1;
        %============================= Plot Grid and Current Light Beam
        if(plotSpec)
            finalfig=figure('name',sprintf('XRF at beam %d with angle %d',i,thetan(n)));
            subplot(1,2,1);
            plotGrid(xc,omega,[m(2) m(1)]);
            hold on;
            plot(DetKnot(:,1),DetKnot(:,2),'k+-',SourceKnot(:,1),SourceKnot(:,2),'m+-',SSDknot(:,1),SSDknot(:,2),'g+-','LineWidth',0.5)
            axis equal 
            W_temp=W; W_temp(:,:,1)=1e-3*W(:,:,1);
            imagesc(xc(:,1),xc(:,2),sum(W_temp,3))
            axis xy
            set(gcf,'Units','normalized')
            set(gca,'Units','normalized')
            ax = axis;
            ap = get(gca,'Position');
            xp = ([SourceKnot(i,1),DetKnot(i,1)]-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
            yp = ([SourceKnot(i,2),DetKnot(i,2)]-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
            ah=annotation('arrow',xp,yp,'Color','r','LineStyle','--');
        end
        %=================================================================
        [index,Lvec,linearInd]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta);
        if(isempty(index))
            EmptyBeam=[EmptyBeam,[n;i]];
        else
            RealBeam=[RealBeam,[n;i]];
        end
        %%%%%%%%================================================================
        GlobalInd{n,i}=index;
        if(synthetic)
            xrfSub=zeros(1,numChannel);
        else
            xrfSub_raw=zeros(1,numChannel_raw);
            xrfSub_decom=zeros(1,numChannel_decom);
        end
        for j=1:size(index,1)
            CurrentCellCenter=[(index(j,1)-1/2)*dz(1)-abs(omega(1)),(index(j,2)-1/2)*dz(2)-abs(omega(3))];
            currentInd=sub2ind(m,index(j,2),index(j,1));
            L(n,i,index(j,2),index(j,1))=Lvec(j);
            % L(sub2ind([numThetan,nTau+1,m],n,i,index(j,2),index(j,1)))=Lvec(j);
            if(j==1)
                I_incident=1;
                temp_sum=0;
                SelfInd{n,i,currentInd}{5}=sub2ind(m,index(end:-1:j+1,2),index(end:-1:j+1,1));
                SelfInd{n,i,currentInd}{6}=kron(MU_e(:,1,1)',Lvec(j));%Lvec(end:-1:j+1));
            elseif(j>1 & j<size(index,1))
                temp_sum=temp_sum+Lvec(j-1)*MU_XTM(index(j-1,2),index(j-1,1));
                I_incident=exp(-temp_sum);
                SelfInd{n,i,currentInd}{1}=sub2ind(m,index(j-1:-1:1,2),index(j-1:-1:1,1));
                SelfInd{n,i,currentInd}{3}=kron(Lvec(j-1:-1:1),MU_e(:,1,1)');
                SelfInd{n,i,currentInd}{5}=sub2ind(m,index(end:-1:j+1,2),index(end:-1:j+1,1));
                SelfInd{n,i,currentInd}{6}=kron(MU_e(:,1,1)',Lvec(j));%Lvec(end:-1:j+1));
            else
                temp_sum=temp_sum+Lvec(j-1)*MU_XTM(index(j-1,2),index(j-1,1));
                I_incident=exp(-temp_sum);
                SelfInd{n,i,currentInd}{1}=sub2ind(m,index(j-1:-1:1,2),index(j-1:-1:1,1));
                SelfInd{n,i,currentInd}{3}=kron(Lvec(j-1:-1:1),MU_e(:,1,1)');
            end
            %% ===========================================================
            Wsub=reshape(W(index(j,2),index(j,1),:),[NumElement,1]);
            %% Self-absorption
            in_after=find(inpolygon(xc(:,1),xc(:,2),[CurrentCellCenter(1) SSDknot(1,1) SSDknot(NumSSDlet,1)],[CurrentCellCenter(2) SSDknot(1,2) SSDknot(NumSSDlet,2)])); 
            in_after=setdiff(in_after,currentInd); %% energy is not attenuated from the source point
            % figure(19);plot(omega([1 2 2 1 1]),omega([3 3 4 4 3]),'r.-');
            % hold on; fill([CurrentCellCenter(1) SSDknot(1,1) SSDknot(NumSSDlet,1)],[CurrentCellCenter(2) SSDknot(1,2) SSDknot(NumSSDlet,2)],'b.-');
            % hold off;
            area_xrf(n,i,currentInd)=prod(dz)*length(in_after);
            SelfInd{n,i,currentInd}{2}=in_after;
            I_after = exp(-sum(MU_after(in_after,:)',2)./(length(in_after)+1)*area_xrf(n,i,currentInd));%% Attenuation of Flourescent energy emitted from current pixel
            %% ====================================================================
            if(synthetic)
                xrfSub=xrfSub+Lvec(j)*I_incident*(I_after.*Wsub)'*M;% fluorescence emitted from current pixel
            else
                xrfSub_raw=xrfSub_raw+Lvec(j)*I_incident*(I_after.*Wsub)'*M_raw;% fluorescence emitted from current pixel
                xrfSub_decom=xrfSub_decom+Lvec(j)*I_incident*(I_after.*Wsub)'*M_decom;% fluorescence emitted from current pixel
           end
        end
        Rdis(i)=I0*exp(-eX'*(MU_XTM.*reshape(L(sub2ind([numThetan,nTau+1,prod(m)],n*ones(1,prod(m)),i*ones(1,prod(m)),1:prod(m))),m))*eY); %% Discrete case
        if(synthetic)
            XRF(n,i,:)=xrfSub;%+0.1*rand(size(xrfSub));
        else
            XRF_raw(n,i,:)=xrfSub_raw;
            XRF_decom(n,i,:)=xrfSub_decom;
        end
        clear xrfSub xrfSub_raw xrfSub_decom
        if(plotSpec)
            figure(finalfig)
            subplot(1,2,2);
            plot(DetChannel,squeeze(XRF_raw(n,i,:)),'r-')
            xlabel('Energy Channel','fontsize',12); ylabel('Intensity','fontsize',12)
            pause(1);
        end
    end
    % L(n,:,:,:)=L(n,:,:,:)./(repmat(sum(L(n,:,:,:),2),[1,nTau+1,1,1])+1)*prod(dz);
    DisR(:,n)=Rdis';
    clear Rdis temp
end
L=L(:);

clear Wsub A knot DetKnot DetKnot0 SSDknot SSDlet SourceKnot SourceKnot0 MU_after
clear xbox ybox eX eY  x y xc 
if(~synthetic)
    XRF_Simulated_raw=XRF_raw;
    XRF_Simulated_decom=XRF_decom;
    DisR_Simulated=DisR;
    XRF_decom=permute(data_xrf_decom(:,:,:),[2 3 1]);
    if(truncChannel)
        % if(strcmp(sample,'Rod'))
        %     truncInd=find(DetChannel> 5 & DetChannel <9);
        % elseif(strcmp(sample,'Seed'))
        %     truncInd=find(DetChannel> 5 & DetChannel <9);
        % end
        DetChannel_raw=DetChannel(truncInd);
        numChannel_raw=length(truncInd);
        XRF_Simulated_raw=XRF_raw(:,:,truncInd);
        XRF_raw=permute(data_xrf_raw(truncInd,:,:),[2 3 1]);
        M_raw=M_raw(:,truncInd);
        clear truncInd DetChannel
    else
        XRF_raw=permute(data_xrf_raw,[2 3 1]);
    end
    DisR=squeeze(sum(data_xrt(:,:,:),1))';
    % save(['Simulated_',num2str(N(1)),'_',num2str(numThetan),'_',num2str(NumElement),'_',num2str(numChannel_raw), sample,'.mat'],'DisR_Simulated','XRF_raw','XRF_Simulated','DisR');
    clear data_xrt data_xrf_decom data_xrf_raw 
    % clear XRF_Simulated DisR_Simulated
end

SigMa_XRF=1;
SigMa_XTM=1;



