%%%Simulate XRF and XRT of a given object with predifined detector and beam
%%% Travelling of fluorescece photons is approximated as numerical disrectized paths in the range of solid angle
global plotSpecSingle BeforeEmit plotTravel SSDlet
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
XRF=zeros(numThetan,nTau+1,numChannel);
XRF_decom=zeros(numThetan,nTau+1,numChannel_decom);
DisR=zeros(nTau+1,numThetan);
L=zeros(numThetan,nTau+1,m(1),m(2));%zeros(numThetan*(nTau+1)*prod(m),1);
if(synthetic)
    fprintf(1,'====== Fluorescence Detector Resolution is %d\n',numChannel);
else
    fprintf(1,'====== Fluorescence Detector Resolution is %d\n',numChannel_raw);
end
% temp_test=zeros(m(1),m(2),NumElement);
% for i=1:NumElement
%     temp_test(:,:,i)=flipud(reshape(MU_after(:,i),m(1),m(2))');
% end
% temp_test=reshape(temp_test,mtol,NumElement);

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
        XRF(n,i,:) = zeros(numChannel,1);
        BeforeEmit=1;
        %============================= Plot Grid and Current Light Beam
        show_beam=23;%floor(nTau/2);
        if(plotSpec & i==show_beam)
            finalfig=figure('name',sprintf('XRF at beam %d with angle %d',i,thetan(n)));
            subplot(1,2,1);
            plotGrid(xc,omega,[m(2) m(1)]);
            hold on;
            plot(DetKnot(:,1),DetKnot(:,2),'k+-',SourceKnot(:,1),SourceKnot(:,2),'m+-',SSDknot(:,1),SSDknot(:,2),'g+-','LineWidth',0.5)
            axis equal 
            imagesc((x(1:end-1)+x(2:end))/2,(y(1:end-1)+y(2:end))/2,sum(W,3))
            plot([SourceKnot(i,1) DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)],'r--')
            axis xy
            % set(gcf,'Units','normalized')
            % set(gca,'Units','normalized')
            % ax = axis;
            % ap = get(gca,'Position');
            % xp = ([SourceKnot(i,1),DetKnot(i,1)]-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
            % yp = ([SourceKnot(i,2),DetKnot(i,2)]-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
            % ah=annotation('arrow',xp,yp,'Color','r','LineStyle','--');
        end
        %=================================================================
        [index,Lvec,linearInd]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta);
        %%%%%%%%================================================================
        xrfSub=zeros(1,numChannel);
        xrfSub_decom=zeros(1,numChannel_decom);
        for j=1:size(index,1)
            CurrentCellCenter=[(index(j,1)-1/2)*dz(1)-abs(omega(1)),(index(j,2)-1/2)*dz(2)-abs(omega(3))];
            currentInd=sub2ind(m,index(j,2),index(j,1));
            L(n,i,index(j,2),index(j,1))=Lvec(j);
            if(j==1)
                I_incident=1;
                temp_sum=0;
            elseif(j>1 & j<size(index,1))
                temp_sum=temp_sum+Lvec(j-1)*MU_XTM(index(j-1,2),index(j-1,1));
                I_incident=exp(-temp_sum);
            else
                temp_sum=temp_sum+Lvec(j-1)*MU_XTM(index(j-1,2),index(j-1,1));
                I_incident=exp(-temp_sum);
            end
            %% ===========================================================
            Wsub=reshape(W(index(j,2),index(j,1),:),[NumElement,1]);
            %% Self-absorption
            if(NoSelfAbsorption)
                I_after=1;
            else
                in_after=find(inpolygon(xc(:,1),xc(:,2),[CurrentCellCenter(1) SSDknot(1,1) SSDknot(NumSSDlet,1)],[CurrentCellCenter(2) SSDknot(1,2) SSDknot(NumSSDlet,2)])); 
                in_after=setdiff(in_after,currentInd); %% energy is not attenuated from the source point
                I_after = exp(-sum(MU_after(in_after,:)',2)./(length(in_after)+1)*prod(dz)*length(in_after));%% Attenuation of Flourescent energy emitted from current pixel
            end
            %% ====================================================================
            xrfSub=xrfSub+Lvec(j)*I_incident*(I_after.*Wsub)'*M_raw;% fluorescence emitted from current pixel
            xrfSub_decom=xrfSub_decom+Lvec(j)*I_incident*(I_after.*Wsub)'*M_decom;% fluorescence emitted from current pixel
        end
        Rdis(i)=I0*exp(-eX'*(MU_XTM.*reshape(L(sub2ind([numThetan,nTau+1,prod(m)],n*ones(1,prod(m)),i*ones(1,prod(m)),1:prod(m))),m))*eY); %% Discrete case
        XRF(n,i,:)=xrfSub;%+0.1*rand(size(xrfSub));
        XRF_decom(n,i,:)=xrfSub_decom;%+0.1*rand(size(xrfSub));
        clear xrfSub
        if(plotSpec & i==show_beam)
            figure(finalfig)
            subplot(1,2,2);
            semilogy(DetChannel,squeeze(XRF(n,i,:)),'r-')
            xlabel('Energy Channel','fontsize',12); ylabel('Intensity','fontsize',12)
            drawnow;
        end
    end
    DisR(:,n)=Rdis';
    clear Rdis temp
end

