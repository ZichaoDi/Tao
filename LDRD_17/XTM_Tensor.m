%%%Simulate XRT of a given object with predifined detector and beam
global plotSpecSingle BeforeEmit N 
global LogScale NonEmptyBeam synthetic 
plotSpecSingle=0;
more off;
%%-------------------------- Set up center of rotation
Tol=1e-2; 
omega=[-2     2    -2     2].*Tol;
m=[current_n current_n]; %Numerical Resolution
dz=[(omega(2)-omega(1))/m(2) (omega(4)-omega(3))/m(1)];

cor=dz(1)/1*(m(1)/3-N(1)/2);
cr=[0 0; -cor -cor; -cor cor; cor cor; cor -cor;cor 0;0 cor;-cor 0;0 -cor];
cr=cr([1 4 6],:);
if(~exist('ind_cr','var'))
  ind_cr=1;
end
Delta_D=[0,1*cor];
delta_d0=Delta_D(1);
%%-----------------------------------------------
Define_Detector_Beam_Gaussian; %% provide the beam source and Detectorlet
DefineObject_Gaussian; % Produce W, MU_XTM
%%==============================================================
Energy=BindingEnergy;
mtol=prod(m);
eX=ones(m(1),1);
eY=ones(m(2),1);
NonEmptyBeam=[];
L=sparse(numThetan*(nTau+1),prod(m));
DisR_Simulated=zeros(nTau+1,numThetan);
fprintf(1,'====== Fluorescence Detector Resolution is %d\n',numChannel);
plotDisBeam=0;
plotRotation=0;
if(plotDisBeam | plotRotation)
    figure,
end
n_delta=2*numThetan;
if(n_delta==2*numThetan)
    rng('default');
    pert=(-1+2*rand(numThetan,2))*dz(1)*2;
    pert=cumsum(pert,1)*1;
elseif(n_delta==4)
    pert1=1*dz(1)*10; pert2=0.5*dz(1)*10;
    pert=[pert1*ones(floor(numThetan/2),2);pert2*ones(numThetan-floor(numThetan/2),2)];
else
    pert=zeros(numThetan,2)*dz(1)*1;
end
xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
if(ind_cr==2)
    delta0=repmat(cr(ind_cr,:),numThetan,1)+pert*1;
    % delta0(:,1)=0;
    % r=dz(1)*5;
    % delta0=[sin(theta/2)*r,cos(theta/2)*r];
    % delta0=[cos(-theta)*r,sin(-theta)*r];
    delta0(delta0==inf)=0;
    % load integerShift;
    % delta0(:,2)=dd1;
    % delta0(:,1)=dd2;
else
    delta0=sparse(numThetan,2);
end
theta=thetan'/180*pi;
% MU_XTM=xs;%imrotate(MU_XTM,60,'crop');
for n=1:numThetan
    if(mod(n,10)==0)
        fprintf(1,'====== Angle Number  %d of %d: %d\n',n,numThetan,thetan(n));
    end
    TransMatrix=[cos(theta(n)) sin(theta(n)) delta0(n,1)-cos(theta(n))*delta0(n,1)-sin(theta(n))*delta0(n,2); -sin(theta(n)) cos(theta(n)) delta0(n,2)+sin(theta(n))*delta0(n,1)-cos(theta(n))*delta0(n,2);  0 0 1];
    DetKnot=TransMatrix*[DetKnot0'; ones(1,size(DetKnot0,1))];
    DetKnot=DetKnot(1:2,:)';
    SourceKnot=TransMatrix*[SourceKnot0'; ones(1,size(DetKnot0,1))];
    SourceKnot=SourceKnot(1:2,:)';
    %% =========================================
    Rdis_true=mean(I0(:))'.*ones(nTau+1,1);
    if(plotDisBeam)
        hold on;
        col=2;
        subplot(numThetan,col,col*(n-1)+1)
        [px,py]=meshgrid(linspace(omega(1),omega(2),m(1)+1),linspace(omega(3),omega(4),m(2)+1));
        imagesc(linspace(omega(1),omega(2),m(1)+1),linspace(omega(3),omega(4),m(2)+1),MU_XTM);
        title(sprintf([num2str(int16(thetan(n))),'%c'], char(176)));
        hold on;
        plot(SourceKnot(:,1),SourceKnot(:,2),'r.-',DetKnot(:,1),DetKnot(:,2),'k.-');
        axis([-0.07 0.07 -0.07 0.07])
        axis square;
        drawnow;
    end
    if(plotRotation)
        output(:,:,:,n)=rotateAround(padarray(W,[ceil(N(1)/2),ceil(N(1)/2)]),delta0(n,1)/dz(1)+N(1),delta0(n,2)/dz(1)+N(1),thetan(n));
        % subplot(5,ceil(numThetan/5),n)
        imagesc(sum(output(:,:,:,n),3)); axis xy image; hold on; plot(delta0(n,1)/dz(1)+N(1),delta0(n,2)/dz(1)+N(1),'r*');
        hold on;
        set(gca,'xtick',[],'ytick',[]);
        drawnow;
    end
    for i=1:nTau+1 %%%%%%%%%========================================================
        % Initialize
        BeforeEmit=1;
        %============================= Plot Grid and Current Light Beam
        %=================================================================
        [index,Lvec,linearInd]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta(n));
        %%%%%%%%================================================================
        if(~isempty(index)& norm(Lvec)>0)
            NonEmptyBeam=[NonEmptyBeam,sub2ind([numThetan,nTau+1],n,i)];
            currentInd=sub2ind(m,index(:,2),index(:,1));
            L(sub2ind([numThetan,nTau+1],n,i),currentInd)=Lvec/sum(Lvec);
            Rdis_true(i)=mean(I0(:))*exp(-eX'*(MU_XTM.*reshape(L(sub2ind([numThetan,nTau+1],n,i),:),m))*eY);%% Discrete case
        end
    end
        DisR_Simulated(:,n)=Rdis_true';
        if(plotDisBeam)
           subplot(numThetan,col,col*(n-1)+2)
           plot(1:nTau+1,DisR_Simulated(:,n),'r.-')
           if(n==1)
               % xlabel('Beamlet','fontsize',12); ylabel('Intensity','fontsize',12)
           end
           drawnow;
        end
end
% scale=1e2;%/sqrt(N(1));
% L=map1D(L*scale,[0,1]);
%%==============================================================
