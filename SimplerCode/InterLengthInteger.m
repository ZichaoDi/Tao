function [DisR,I,Ltol,thetan,omega,m,dTau]=InterLength;%(xc,theta,dTau)
%%% L: the matrix of intersection length of beam line travelling
%%% throug the object at angle theta
%%% xc: discrete grids on O
%%% dTau*sin(theta): width of each detectorlet
%%% the beam line is y=x*tan(theta)+t*dTau, t=1,...,N_d
close all
%%==== Discrete Test object
% I=zeros(100,100);
% I(25:75,25:75)=1;
% % I=imrotate(I,30);
% % I=double(imread('LoganReference.tiff'));
%%%%%%%%%%%%%%%%%%%%=============================================
Dis=10;
x=linspace(0,Dis,Dis);
y=linspace(0,Dis,Dis);
[X,Y] = meshgrid(x,y);
center=size(X)./2;
I=(X-center(1)).^2+(Y-center(2)).^2;
figure('name','Object');
% I=GrayScale(I);
imagesc(I); axis xy
e=ones(size(I,1),1);
m=size(I);
omega=[0 Dis 0 Dis];

xc = getNodalGrid(omega,m);
figure('name','grids & beam line')
plotGrid(xc,omega,m);
%%==== Intersection
dTau=1;
thetan=1:2:180;
DisR=[];
ConR=[];
RadonR=[];
nTau=ceil(sqrt(omega(2)^2+omega(4)^2)/dTau);
Ltol=cell(length(thetan),nTau);
for n=68;%1:length(thetan)
    theta=thetan(n)/180*pi;
    
    z=(0:omega(2))';
    if(theta<=pi/2)
    diagL1=ceil(omega(2)*tan(theta)/dTau);
    diagL2=ceil(omega(4)/dTau);
    else
    diagL1=0;%ceil(omega(2)*tan(theta)/dTau);
    diagL2=ceil((omega(4)-omega(2)*tan(theta))/dTau);
    end
    Ttol=linspace(-diagL1+0,diagL2-0,nTau);
    for i=1:length(Ttol)
        t=Ttol(i);
            Q=[[z, tan(theta).*z+t*dTau];[(z-t*dTau)./tan(theta),z]];

        indx=find(Q(:,1)<omega(1) |Q(:,1)>omega(2)); indy=find(Q(:,2)<omega(3) |Q(:,2)>omega(4));
        Q=setdiff(Q,Q([indx;indy],:),'rows');
        Q=unique(Q,'rows');
        Lvec=sqrt((Q(2:end,1)-Q(1:end-1,1)).^2+(Q(2:end,2)-Q(1:end-1,2)).^2);
        if(theta<=pi/2)
        index=floor(Q+1);
        else
            index=[ceil(Q(:,2)),floor(Q(:,1))+1];
        end
        inx=find(index(:,1)<=omega(1) |index(:,1)>m(1)); iny=find(index(:,2)<=omega(3) |index(:,2)>m(2));
        index=setdiff(index,index([inx,iny],:),'rows');
        %%%%%%%%%================================================================
        if(theta<=pi/2)
        Ax1=[-t*dTau/tan(theta), omega(2)];
        Ay1=[0 tan(theta).*omega(2)+t*dTau];
        else
            Ax1=[-t*dTau/tan(theta), omega(1)];
        Ay1=[0 tan(theta).*omega(1)+t*dTau];
        end
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        [Ax, Ay] = polyxpoly(Ax1, Ay1, xbox, ybox);
        %         Ax=unique(sort(Ax));Ay=unique(sort(Ay));
        [Ax,in]=sort(Ax);A=unique([Ax,Ay(in)],'rows');Ax=A(:,1);Ay=A(:,2);
        if(isempty(Ax) | (length(Ax)==1 &length(Ay)==1 ))
%             fprintf('no intersection \n')
            L=sparse(zeros(m(1),m(2)));
            Rdis(i)=0;
            R(i)=0;
        else
        drawnow;
        hold on; plot(Ax1,Ay1,'r.-');
        if(theta<=pi/2)
                    plot(index(:,1)-0.5,index(:,2)-0.5,'bo',Q(:,1),Q(:,2),'g.');
        else
        plot((index(:,2)-0.5),(index(:,1)-0.5),'bo',Q(:,1),Q(:,2),'g.');
        end
        axis(omega);
        pause;
        %%%%%%%%%================================================================
        L=zeros(m(1),m(2));
        for j=1:size(index,1)
            L(index(j,1),index(j,2))=Lvec(j);
        end
        Ltol{n,i}=L;
        Rdis(i)=e'*(I.*L)*e; %% Discrete case
        %%============ Continuous integral
        syms u;
        R(i)=int(((Ax(1)+(Ax(2)-Ax(1))*u-center(1))^2+(Ay(1)+(Ay(2)-Ay(1))*u-center(2))^2)*sqrt((Ax(2)-Ax(1))^2+(Ay(2)-Ay(1))^2),0,1);
        %%==========================================================
        end
    end
    DisR(:,n)=Rdis';
    ConR(:,n)=R';
    RadonR(:,n)=radon(I,theta/pi*180);
end
iDisR=iradon(DisR,thetan);
iConR=iradon(ConR,thetan);
iRadonR=iradon(RadonR,thetan);
% figure('name','Projection at angle theta')
% Rspace=linspace(Ttol(1),Ttol(end),length(RadonR));
% plot(Ttol,R,'r.-',Ttol, Rdis,'g*-',Rspace,RadonR,'bo-');
figure('name','Discrete Reconstruction');imagesc(iDisR);
figure('name','Continuous Reconstruction');imagesc(iConR);
figure('name','Matlab Reconstruction');imagesc(iRadonR);

