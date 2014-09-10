 function [W,XRF,RMlocal,I,M,Energy,EnergyChannel,Ltol,GlobalInd,thetan,m,nTau]=XRFreconCompond;
%%% L: the matrix of intersection length of beam line travelling
%%% throug the object at angle theta
%%% xc: discrete grids on O
%%% dTau: width of each detectorlet
%%% the beam line is y=x*tan(theta)+t*dTau/cos(theta), t=1,...,N_d
global DetChannel
more off;
close all
%%==== Discrete Test object
% I=zeros(100,100);
% I(25:75,25:75)=1;
% % I=imrotate(I,30);
% % I=double(imread('LoganReference.tiff'));
%%%%%%%%%%%%%%%%%%%%=============================================


% Do you want to see the spectra? If so plotflag = 1
plotflag = 0;

% State whether you are on Wendy's computer:
wendy = 0;
Dis=4;
omega=[-Dis/2 Dis/2 -Dis/2 Dis/2];
x=linspace(omega(1),omega(2),Dis);
y=linspace(omega(3),omega(4),Dis);
[X,Y] = meshgrid((x(1:end-1)+x(2:end))./2,(y(1:end-1)+y(2:end))./2);
center=[0 0];
 xrfIntensity_mat;
m=size(I);
% [xrf_intensity,beam_energy]=xrfIntensity_mat;
I0=beam_energy;
M=xrf_intensity(:,2);
Energy=xrf_intensity(:,1);
EnergyChannel=xrf_intensity(:,1);
NumElement=size(M,1);
W=zeros(m(1),m(2),NumElement); %weight fraction of pixel for each element
O=zeros(m(1),m(2),NumElement); % O=W.*Z for plotting
Z=[30 20 35 26]; % rename as AtomicNum
%%% fill comments later
W(1,1,:)=[1 0 0 0];
W(1,2,:)=[0 1 0 0];
W(1,3,:)=[1 0 0 0];
W(2,1,:)=[0 1 0 0];
W(2,2,:)=   [0 0 0.5 0.5] ;
W(2,3,:)=  [0 1 0 0];
W(3,1,:)=   [0 0 1 0] ;
W(3,2,:)=   [1 0 0 0] ;
W(3,3,:)=  [0 1 0 0]; 
W(4,4,1) =1;
W(1:3,4,2) =1;
W(4,1:3,3) =1;
O=zeros(4,4,NumElement); % O=W.*Z for plotting

%%% fill comments l
for e=1:4
    O(:,:,e)=W(:,:,e).*Z(e);
end
figure('name','AttenuationMap');
imagesc(I); axis xy
figure('name','Object');
for e=1:4
    subplot(2,2,e);
    image(O(:,:,e));
    axis xy
end


dz=(Dis)/m(1);
alpha=atan(m(2)/m(1));
%%%================================================================
xc = getNodalGrid(omega,m);
fig=figure('name','grids & beam line');
plotGrid(xc,omega,m); hold on;
fig1=[]; fig2=[]; fig3=[];
%%==== Intersection
dTau=dz;
thetan=1:2:180;%[1  45 91];%
Tau=sqrt(m(1)^2+m(2)^2);
nTau=ceil(Tau/dTau)+2; % number of detectorlets
DisNtau=nTau;
tol=10;
detS0=[Tau/2*tan(alpha)+tol, ceil(Tau/2)]; %initiate detector location
detE0=[Tau/2*tan(alpha)+tol,-ceil(Tau/2)];
knot=linspace(detS0(2),detE0(2),nTau+1)';%% detectorlet knot points
DetKnot0=[repmat(detS0(1),size(knot)),knot];
SourceS0=[-Tau/2*tan(alpha)-tol, ceil(Tau/2)];
SourceE0=[-Tau/2*tan(alpha)-tol,-ceil(Tau/2)];
knot=linspace(SourceS0(2),SourceE0(2),nTau+1)';%% detectorlet knot points
SourceKnot0=[repmat(SourceS0(1),size(knot)),knot];

XRF=cell(length(thetan),nTau+1);

Ltol=cell(length(thetan),nTau+1);
GlobalInd=cell(length(thetan),nTau+1);
RMlocal=zeros(m(1),m(2),NumElement);
fprintf(1,'====== Detector Resolution is %d\n',nTau);
for n=1:length(thetan)
    
    fprintf(1,'====== Angle Number  %d of %d\n',n,length(thetan))
    
    theta=thetan(n)/180*pi;
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    detS=detS0*TransMatrix;
    detE=detE0*TransMatrix;
    SourceS=SourceS0*TransMatrix;
    SourceE=SourceE0*TransMatrix;
DetKnot=DetKnot0*TransMatrix;
SourceKnot=SourceKnot0*TransMatrix;
    if wendy && plotflag
        figure(fig);
        set(fig3,'visible','off');
        for i=1:nTau+1
            fig3=plot([SourceKnot(i,1),DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)],'r.-',SourceKnot(i,1),SourceKnot(i,2),'bo');hold on;
        end
        pause;
    end
    
    
    for i=1:DisNtau+1 %%%%%%%%%================================================================
        % Initialize
        XRF{n,i} = zeros(numChannel,1);
         
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        if wendy % Change this later, for now we don't have vpa :(
            [Ax, Ay] = polyxpoly([double(vpa(SourceKnot(i,1),5)),DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)], xbox, ybox);
        else
            [Ax, Ay] = polyxpoly([SourceKnot(i,1),DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)], xbox, ybox);
        end
        if(isempty(Ax) | (length(Ax)==1 &length(Ay)==1 ))
            %             fprintf('no intersection \n')
            index=[];
            Q=[];
            L=sparse(zeros(m(1),m(2)));
            
        else
            if(theta<=pi/2)
                [Ax,in]=sort(Ax);A=unique([Ax,Ay(in)],'rows');Ax=A(:,1);Ay=A(:,2);
            else
                [Ay,in]=sort(Ay);
                Asub=[Ax(in),Ay];
                [A,in]=unique(Asub,'rows');
                if(size(A,1)>1)
                    A=Asub(in,:);Ax=A(:,1);Ay=A(:,2);
                end
            end
            
            if(theta==pi/2)
                Q=[repmat(Ax(1),size(y')),y'];
            elseif(theta==0 | theta==pi)
                Q=[x',repmat(Ay(1),size(x'))];
            else
                %%%%%%%%%================================================================
                Q=[[x', (Ay(2)-Ay(1))/(Ax(2)-Ax(1)).*x'+(Ay(1)*Ax(2)-Ay(2)*Ax(1))/(Ax(2)-Ax(1))];...
                    [(y'-(Ay(1)*Ax(2)-Ay(2)*Ax(1))/(Ax(2)-Ax(1)))./((Ay(2)-Ay(1))/(Ax(2)-Ax(1))),y']];
            end
            indx=find(Q(:,1)<omega(1) |Q(:,1)>omega(2)); indy=find(Q(:,2)<omega(3) |Q(:,2)>omega(4));
            Q=setdiff(Q,Q([indx;indy],:),'rows');
            Q=unique(Q,'rows');
            Lvec=sqrt((Q(2:end,1)-Q(1:end-1,1)).^2+(Q(2:end,2)-Q(1:end-1,2)).^2);
            %         Lvec=Lvec(find(Lvec>1e-10));
            %%%%%%%%%================================================================
            %%%%%%%%%================================================================
            if(theta<=pi/2)
                if wendy % Change this later, for now we don't have vpa :(
                    index=[floor(double(vpa((Q(:,1)+abs(omega(1)))/dz,10)))+1,floor(double(vpa((Q(:,2)+abs(omega(3)))/dz,10)))+1];
                else
                    index=[floor(double((Q(:,1)+abs(omega(1)))/dz))+1,floor(double((Q(:,2)+abs(omega(3)))/dz))+1];
                end
            else
                if wendy % Change this later, for now we don't have vpa :(
                    index=[floor(double(vpa((Q(:,1)+abs(omega(1)))/dz,10)))+1,ceil(double(vpa((Q(:,2)+abs(omega(3)))/dz,10)))];
                else
                    index=[floor(double((Q(:,1)+abs(omega(1)))/dz))+1,ceil(double((Q(:,2)+abs(omega(3)))/dz))];
                end
            end
            index=index(find(index(:,1)>0 & index(:,1)<=m(1)& index(:,2)<=m(1) & index(:,2)>0),:);
            GlobalInd{n,i}=index;
            %%%%%%%%%================================================================
            if plotflag
                figure(fig);
                set(fig1,'visible','off');
                set(fig2,'visible','off');
                drawnow;
                fig1=plot(Ax,Ay,'r.-');%,SourceKnot(i,1),SourceKnot(i,2),'ro',DetKnot(i,1),DetKnot(i,2),'r*');
                if(~isempty(index))
                    %                 if(theta==0 |theta==pi |theta==pi/2)
                    %                     fig2=plot((index(:,2)-0.5)*dz,(index(:,1)-0.5)*dz,'bo');
                    %
                    %                 else
                    fig2=plot((index(:,1)-abs(omega(1)))*dz,(index(:,2)-abs(omega(3)))*dz,'bo',Q(:,1),Q(:,2),'g.');
                    %                 end
                end
            end
            %%%%%%%%%================================================================
            
            L=sparse(zeros(m(1),m(2)));
            RM=cell(m(1),m(2));
            xrfSub=zeros(size(EnergyChannel));
            for j=1:size(index,1)
                L(index(j,2),index(j,1))=Lvec(j);
                if(j==1)
                    I_incident=1;
                    temp_sum=0;
                else
                    temp_sum=temp_sum+L(index(j-1,2),index(j-1,1))*I(index(j-1,2),index(j-1,1));
                    I_incident=exp(-temp_sum);
                end
                Wsub=reshape(W(index(j,2),index(j,1),:),size(M));
                RM{index(j,2),index(j,1)}=Lvec(j)*I_incident*(Wsub.*M);
                temp=zeros(size(RMlocal(index(j,2),index(j,1),:)));
                temp(1,1,:)=Lvec(j)*I_incident*(Wsub.*M);
                RMlocal(index(j,2),index(j,1),:)=RMlocal(index(j,2),index(j,1),:)+temp;
                ExistInd=find(Wsub~=0);
                for tsub=1:length(ExistInd)
                    DetInd=find(Energy(ExistInd(tsub))==EnergyChannel);
                    
                    xrfSub(DetInd)=xrfSub(DetInd)+RM{index(j,2),index(j,1)}(ExistInd(tsub));
                end
            end
            Ltol{n,i}=L;            
            
            p1=GaussianFit1(EnergyChannel,xrfSub);
            if(plotflag)  
                finalfig=figure('name',sprintf('XRF at beam %d with angle %d',i,thetan(n)));  
                plot(DetChannel,p1,'r-')
                xlabel('Energy Channel','fontsize',12); ylabel('Intensity','fontsize',12)
                 pause;
            end
            XRF{n,i}=p1;
        end
        
    end
end

