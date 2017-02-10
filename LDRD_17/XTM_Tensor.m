%%%Simulate XRT of a given object with predifined detector and beam
global plotSpecSingle BeforeEmit N 
global LogScale EmptyBeam synthetic 
plotSpecSingle=0;
more off;
%%-------------------------- Set up center of rotation
Tol=1e-2; 
omega=[-2     2    -2     2].*Tol;
m=[current_n current_n]; %Numerical Resolution
dz=[(omega(2)-omega(1))/m(2) (omega(4)-omega(3))/m(1)];

cor=dz(1)/1*(m(1)/3-N/2);
cr=[0 0; -cor -cor; -cor cor; cor cor; cor -cor;cor 0;0 cor;-cor 0;0 -cor];
cr=cr([1 4 6],:);
if(~exist('ind_cr','var'))
  ind_cr=1;
end
delta0(1)=cr(ind_cr,1); 
delta0(2)=cr(ind_cr,2);
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
EmptyBeam=[];
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
    pert=(-1+2*rand(numThetan,2))*dz(1)*5;
    pert=cumsum(pert,1);
elseif(n_delta==4)
    pert1=1*dz(1)*10; pert2=0.5*dz(1)*10;
    pert=[pert1*ones(floor(numThetan/2),2);pert2*ones(numThetan-floor(numThetan/2),2)];
else
    pert=zeros(numThetan,2);
end
for n=1:numThetan
    if(mod(n,10)==0)
        fprintf(1,'====== Angle Number  %d of %d: %d\n',n,numThetan,thetan(n));
    end
    theta = thetan(n)/180*pi;
    if(ind_cr==2)
        delta0=cr(ind_cr,:)+pert(n,:);
    end
    TransMatrix=[cos(theta) sin(theta) delta0(1)-cos(theta)*delta0(1)-sin(theta)*delta0(2); -sin(theta) cos(theta) delta0(2)+sin(theta)*delta0(1)-cos(theta)*delta0(2);  0 0 1];
    DetKnot=TransMatrix*[DetKnot0'; ones(1,size(DetKnot0,1))];
    DetKnot=DetKnot(1:2,:)';
    SourceKnot=TransMatrix*[SourceKnot0'; ones(1,size(DetKnot0,1))];
    SourceKnot=SourceKnot(1:2,:)';
    %% =========================================
    Rdis_true=mean(I0(:))'.*ones(nTau+1,1);
    xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
    ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
    if(plotDisBeam)
        hold on;
        col=2;
        subplot(numThetan,col,col*(n-1)+1)
        [px,py]=meshgrid(linspace(omega(1),omega(2),m(1)+1),linspace(omega(3),omega(4),m(2)+1));
        imagesc(linspace(omega(1),omega(2),m(1)+1),linspace(omega(3),omega(4),m(2)+1),MU_XTM);
        title(['Angle',num2str(thetan(n))]);
        hold on;
        plot(SourceKnot(:,1),SourceKnot(:,2),'r.-',DetKnot(:,1),DetKnot(:,2),'k.-');
        axis([-0.05 0.05 -0.05 0.05])
        drawnow;
    end
    if(plotRotation)
        output=rotateAround(padarray(W,[N/2,N/2]),delta0(1)/dz(1)+N,delta0(2)/dz(1)+N,thetan(n));
        subplot(5,7,n)
        imagesc(sum(output,3)); hold on; plot(delta0(1)/dz(1)+N,delta0(2)/dz(1)+N,'r*');
        set(gca,'xtick',[],'ytick',[]);
        drawnow;
    end
    for i=1:nTau+1 %%%%%%%%%========================================================
        % Initialize
        BeforeEmit=1;
        %============================= Plot Grid and Current Light Beam
        %=================================================================
        %=================================================================
        [index,Lvec,linearInd]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta);
        %%%%%%%%================================================================
        if(~isempty(index)& norm(Lvec)>0)
            EmptyBeam=[EmptyBeam,(n-1)*numThetan+i];
            currentInd=sub2ind(m,index(:,2),index(:,1));
            L(sub2ind([numThetan,nTau+1],n,i),currentInd)=Lvec;
            Rdis_true(i)=mean(I0(:))*exp(-eX'*(MU_XTM.*reshape(L(sub2ind([numThetan,nTau+1],n,i),:),m))*eY);%% Discrete case
        end
    end
        DisR_Simulated(:,n)=Rdis_true';
        if(plotDisBeam)
           subplot(numThetan,col,col*(n-1)+2)
           plot(1:nTau+1,DisR_Simulated(:,n),'r.-')
           xlabel('Beamlet','fontsize',12); ylabel('Intensity','fontsize',12)
           drawnow;
        end
end
%%==============================================================
