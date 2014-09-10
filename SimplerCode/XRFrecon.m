
%%% L: the matrix of intersection length of beam line travelling
%%% throug the object at angle theta
%%% xc: discrete grids on O
%%% dTau: width of each detectorlet
%%% the beam line is y=x*tan(theta)+t*dTau/cos(theta), t=1,...,N_d
clear all
more off;
close all
time=cputime;
%%==== Discrete Test object
% I=zeros(100,100);
% I(25:75,25:75)=1;
% % I=imrotate(I,30);
% % I=double(imread('LoganReference.tiff'));
%%%%%%%%%%%%%%%%%%%%=============================================

Dis=4;
omega=[-Dis/2 Dis/2 -Dis/2 Dis/2];
x=linspace(omega(1),omega(2),Dis);
y=linspace(omega(3),omega(4),Dis);
[X,Y] = meshgrid((x(1:end-1)+x(2:end))./2,(y(1:end-1)+y(2:end))./2);
center=[0 0];
I=(X-center(1)).^2+(Y-center(2)).^2;
[xrf_intensity,beam_energy]=xrfIntensity;
I0=beam_energy;
M=[xrf_intensity(2,2) xrf_intensity(3,2) xrf_intensity(2,2);...
    xrf_intensity(3,2) xrf_intensity(1,2) xrf_intensity(3,2);...
    xrf_intensity(2,2) xrf_intensity(3,2) xrf_intensity(2,2)];
Energy=[xrf_intensity(2,1) xrf_intensity(3,1) xrf_intensity(2,1);...
    xrf_intensity(3,1) xrf_intensity(1,1) xrf_intensity(3,1);...
    xrf_intensity(2,1) xrf_intensity(3,1) xrf_intensity(2,1)];

EnergyChannel=xrf_intensity(:,1);

figure('name','Object');
imagesc(I); axis xy
e=ones(size(I,1),1);
m=size(I);
dz=(Dis)/m(1);
TolCenter=[(omega(2)+omega(1))/2,(omega(3)+omega(4))/2];
alpha=atan(m(2)/m(1));
%%%================================================================
xc = getNodalGrid(omega,m);
fig=figure('name','grids & beam line');
plotGrid(xc,omega,m); hold on;
%%==== Intersection
dTau=dz;
thetan=[45];%;%1:2:180;%[1  45 91];%
Tau=sqrt(m(1)^2+m(2)^2);
nTau=ceil(Tau/dTau)+30; % number of detectorlets
tol=10;
detS0=[Tau/2*tan(alpha)+tol, ceil(Tau/2)]; %initiate detector location
detE0=[Tau/2*tan(alpha)+tol,-ceil(Tau/2)];
SourceS0=[-Tau/2*tan(alpha)-tol, ceil(Tau/2)];
SourceE0=[-Tau/2*tan(alpha)-tol,-ceil(Tau/2)];

XRF=cell(length(thetan),nTau+1);
fig1=[]; fig2=[]; fig3=[];
RMtot=sparse(zeros(m(1),m(2)));
for n=1:length(thetan)
    theta=thetan(n)/180*pi;
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    detS=detS0*TransMatrix;
    detE=detE0*TransMatrix;
    SourceS=SourceS0*TransMatrix;
    SourceE=SourceE0*TransMatrix;
    knot=linspace(detS(1),detE(1),nTau+1)';%% detectorlet knot points
    DetKnot=[knot,knot.*(detE(2)-detS(2))/(detE(1)-detS(1))+(detE(1)*detS(2)-detS(1)*detE(2))/(detE(1)-detS(1))];
    knot=linspace(SourceS(1),SourceE(1),nTau+1)';%% detectorlet knot points
    SourceKnot=[knot,knot.*(SourceE(2)-SourceS(2))/(SourceE(1)-SourceS(1))+(SourceE(1)*SourceS(2)-SourceS(1)*SourceE(2))/(SourceE(1)-SourceS(1))];
    %     figure(fig);
    %     set(fig3,'visible','off');
    %     for i=1:nTau+1
    %      fig3=plot([SourceKnot(i,1),DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)],'r.-',SourceKnot(i,1),SourceKnot(i,2),'bo');hold on;
    %     end
    %  pause;
    for i=1:nTau+1
        %%%%%%%%%================================================================
        xrfSub=zeros(size(EnergyChannel));
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        [Ax, Ay] = polyxpoly([double(vpa(SourceKnot(i,1),5)),DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)], xbox, ybox);
        if(theta<=pi/2)
            [Ax,in]=sort(Ax);A=unique([Ax,Ay(in)],'rows');Ax=A(:,1);Ay=A(:,2);
        else
            [Ay,in]=sort(Ay);[A,in]=unique([Ax(in),Ay],'rows');
            if(size(A,1)>1)
                A=A(in,:);Ax=A(:,1);Ay=A(:,2);
            end
        end
        
        if(isempty(Ax) | (length(Ax)==1 &length(Ay)==1 ))
            fprintf('no intersection \n')
            index=[];
            Q=[];
            L=sparse(zeros(m(1),m(2)));
        else
            %%%%%%%%%================================================================
            Q=[[x', (Ay(2)-Ay(1))/(Ax(2)-Ax(1)).*x'+(Ay(1)*Ax(2)-Ay(2)*Ax(1))/(Ax(2)-Ax(1))];...
                [(y'-(Ay(1)*Ax(2)-Ay(2)*Ax(1))/(Ax(2)-Ax(1)))./((Ay(2)-Ay(1))/(Ax(2)-Ax(1))),y']];
            indx=find(Q(:,1)<omega(1) |Q(:,1)>omega(2)); indy=find(Q(:,2)<omega(3) |Q(:,2)>omega(4));
            Q=setdiff(Q,Q([indx;indy],:),'rows');
            Q=unique(Q,'rows');
            Lvec=sqrt((Q(2:end,1)-Q(1:end-1,1)).^2+(Q(2:end,2)-Q(1:end-1,2)).^2);
            %         Lvec=Lvec(find(Lvec>1e-10));
            %%%%%%%%%================================================================
            %%%%%%%%%================================================================
            if(theta<=pi/2)
                index=[floor(double(vpa((Q(:,1)+abs(omega(1)))/dz,10)))+1,floor(double(vpa((Q(:,2)+abs(omega(3)))/dz,10)))+1];
                
            else
                index=[floor(double(vpa((Q(:,1)+abs(omega(1)))/dz,10)))+1,ceil(double(vpa((Q(:,2)+abs(omega(3)))/dz,10)))];
                
            end
            index=index(find(index(:,1)>0 & index(:,1)<=m(1)& index(:,2)<=m(1) & index(:,2)>0),:);
            %         index=unique(index,'rows');
            %%%%%%%%%================================================================
            figure(fig);
            set(fig1,'visible','off');
            set(fig2,'visible','off');
            drawnow;
            fig1=plot(Ax,Ay,'r.-',SourceKnot(i,1),SourceKnot(i,2),'ro',DetKnot(i,1),DetKnot(i,2),'r*');
            if(~isempty(index))
                if(theta==0 |theta==pi |theta==pi/2)
                    fig2=plot((index(:,2)-0.5)*dz,(index(:,1)-0.5)*dz,'bo');
                    
                else
                    fig2=plot((index(:,1)-abs(omega(1)))*dz,(index(:,2)-abs(omega(3)))*dz,'bo',Q(:,1),Q(:,2),'g.');
                end
                % axis(omega);
            end
            %%%%%%%%%================================================================
            L=sparse(zeros(m(1),m(2)));
            RM=sparse(zeros(m(1),m(2)));
            for j=1:size(index,1)
                L(index(j,2),index(j,1))=Lvec(j);
                if(j==1)
                    I_incident=1;
                    temp_sum=0;
                else
                    
                    temp_sum=temp_sum+L(index(j-1,2),index(j-1,1))*I(index(j-1,2),index(j-1,1));
                    I_incident=exp(-temp_sum);
                end
                RM(index(j,2),index(j,1))=Lvec(j)*I_incident*M(index(j,2),index(j,1));
                RMtot(index(j,2),index(j,1))=RMtot(index(j,2),index(j,1))+RM(index(j,2),index(j,1));
                DetInd=find(Energy(index(j,2),index(j,1))==EnergyChannel);
                xrfSub(DetInd)=xrfSub(DetInd)+RM(index(j,2),index(j,1));
            end
            XRF{n,i}=[EnergyChannel,xrfSub];
            finalfig=figure('name',sprintf('XRF at beam %d with angle %d',i,thetan(n)));
            [x1,p1]=GaussianFit1(EnergyChannel,xrfSub,I0,finalfig);
            pause;
            
        end
        
    end
    
end


