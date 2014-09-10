function [ConR,I,Ltol,thetan,omega,m,dTau]=InterLength;%(xc,theta,dTau)
%%% L: the matrix of intersection length of beam line travelling
%%% throug the object at angle theta
%%% xc: discrete grids on O
%%% dTau*sin(theta): width of each detectorlet
%%% the beam line is y=x*tan(theta)+t*dTau, t=1,...,N_d
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
Dis=20;
x=linspace(0,Dis,Dis+2);
y=linspace(0,Dis,Dis+2);
[X,Y] = meshgrid((x(1:end-1)+x(2:end))./2,(y(1:end-1)+y(2:end))./2);
center=[Dis/2+0,Dis/2-0];
I=(X-center(1)).^2+(Y-center(2)).^2;
figure('name','Object');
% I=GrayScale(I);
imagesc(I); axis xy
e=ones(size(I,1),1);
m=size(I);
dz=Dis/m(1);
omega=[0 Dis 0 Dis];

xc = getNodalGrid(omega,m);
fig=figure('name','grids & beam line');
plotGrid(xc,omega,m); hold on;
%%==== Intersection
dTau=1;
thetan=1:2:180;%[1 30 45 100];%
DisR=[];
ConR=[];
RadonR=[];
nTau=ceil(sqrt(omega(2)^2+omega(4)^2)/dTau);
Ltol=cell(length(thetan),nTau);
fig1=[]; fig2=[];
for n=1:length(thetan)
    theta=thetan(n)/180*pi;
    
    if(theta==0 | theta==pi)
        diagL1=ceil(omega(1));
        diagL2=ceil(omega(2));
    elseif(theta==pi/2)
        diagL1=ceil(omega(3));
        diagL2=ceil(omega(4));
    elseif(theta<pi/2)
        diagL1=ceil(omega(2)*tan(theta)/dTau);
        diagL2=ceil(omega(4)/dTau);
    else
        diagL1=0;%ceil(omega(2)*tan(theta)/dTau);
        diagL2=ceil((omega(4)-omega(2)*tan(theta))/dTau);
    end
    %     fprintf('at angle %d, the length of the detector is %d \n',thetan(n),diagL2+diagL1);
    Ttol=linspace(-diagL1+0,diagL2-0,nTau);
    for i=1:length(Ttol)
        %%%%%%%%%================================================================
        
        %%%%%%%%%================================================================
        t=Ttol(i);
        if(theta==0 | theta==pi)
            Ax=omega(1:2)';Ay=[t*dTau;t*dTau];
            index=[floor(t*dTau/dz)*ones(m(1),1)+1,[1:m(2)]'];
            inx=find(index(:,1)<=omega(1) |index(:,1)>m(1)); iny=find(index(:,2)<=omega(3) |index(:,2)>m(2));
            index=setdiff(index,index([inx,iny],:),'rows');
            L=sparse(zeros(m(1),m(2)));
            for j=1:size(index,1)
                L(index(j,1),index(j,2))=dz;
            end
            
            Ltol{n,i}=L;
            Rdis(i)=e'*(I.*L)*e; %% Discrete case
            %%============ Continuous integral
%             syms u;
%             R(i)=double(int(((Ax(1)+(Ax(2)-Ax(1))*u-center(1))^2+(Ay(1)+(Ay(2)-Ay(1))*u-center(2))^2)*sqrt((Ax(2)-Ax(1))^2+(Ay(2)-Ay(1))^2),0,1));
                    R(i)=integral(@(u)((Ax(1)+(Ax(2)-Ax(1)).*u-center(1)).^2+(Ay(1)+(Ay(2)-Ay(1)).*u-center(2)).^2).*sqrt((Ax(2)-Ax(1)).^2+(Ay(2)-Ay(1)).^2),0,1);
        elseif(theta==pi/2)
            Ax=[t*dTau;t*dTau];Ay=omega(3:4)';
            index=[[1:m(1)]',ceil(t*dTau/dz)*ones(m(2),1)];
            inx=find(index(:,1)<=omega(1) |index(:,1)>m(1)); iny=find(index(:,2)<=omega(3) |index(:,2)>m(2));
            index=setdiff(index,index([inx,iny],:),'rows');
            L=sparse(zeros(m(1),m(2)));
            for j=1:size(index,1)
                L(index(j,2),index(j,1))=dz;
            end
            
            Ltol{n,i}=L;
            Rdis(i)=e'*(I.*L)*e; %% Discrete case
            %%============ Continuous integral
%             syms u;
%             R(i)=double(int(((Ax(1)+(Ax(2)-Ax(1))*u-center(1))^2+(Ay(1)+(Ay(2)-Ay(1))*u-center(2))^2)*sqrt((Ax(2)-Ax(1))^2+(Ay(2)-Ay(1))^2),0,1));
            
            R(i)=integral(@(u)((Ax(1)+(Ax(2)-Ax(1)).*u-center(1)).^2+(Ay(1)+(Ay(2)-Ay(1)).*u-center(2)).^2).*sqrt((Ax(2)-Ax(1)).^2+(Ay(2)-Ay(1)).^2),0,1);      
        else
            Q=[[x', tan(theta).*x'+t*dTau];[(y'-t*dTau)./tan(theta),y']];
            indx=find(Q(:,1)<omega(1) |Q(:,1)>omega(2)); indy=find(Q(:,2)<omega(3) |Q(:,2)>omega(4));
            Q=setdiff(Q,Q([indx;indy],:),'rows');
            Q=unique(Q,'rows');
            Lvec=sqrt((Q(2:end,1)-Q(1:end-1,1)).^2+(Q(2:end,2)-Q(1:end-1,2)).^2);
            Lvec=Lvec(find(Lvec>1e-10));
            %%%%%%%%%================================================================
            %%%%%%%%%================================================================
            if(theta<=pi/2)
                index=floor(myvpa(Q/dz+1));
            else
                index=[floor(myvpa(Q(:,1)./dz))+1,ceil(myvpa(Q(:,2)./dz))];
                
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
            if(~isempty(Ax))
                if(theta<=pi/2)
                    [Ax,in]=sort(Ax);
                    A=unique([Ax,Ay(in)],'rows');
                    Ax=A(:,1);Ay=A(:,2);
                else
                    [Ay,in]=sort(Ay);[A,in]=unique([Ax(in),Ay],'rows');
                    if(size(A,1)>1)
                        A=A(in,:);Ax=A(:,1);Ay=A(:,2);
                    end
                end
            end
            
            if(isempty(Ax) | (length(Ax)==1 &length(Ay)==1 ))
                %             fprintf('no intersection \n')
                L=sparse(zeros(m(1),m(2)));
                Rdis(i)=0;
                R(i)=0;
                Ltol{n,i}=L;
            else
                
                L=sparse(zeros(m(1),m(2)));
                for j=1:size(index,1)
                    L(index(j,2),index(j,1))=Lvec(j);
                end
                
                Ltol{n,i}=L;
                Rdis(i)=e'*(I.*L)*e; %% Discrete case
                %%============ Continuous integral
%                 syms u;
%                 R(i)=double(int(((Ax(1)+(Ax(2)-Ax(1))*u-center(1))^2+(Ay(1)+(Ay(2)-Ay(1))*u-center(2))^2)*sqrt((Ax(2)-Ax(1))^2+(Ay(2)-Ay(1))^2),0,1));
          R(i)=integral(@(u)((Ax(1)+(Ax(2)-Ax(1)).*u-center(1)).^2+(Ay(1)+(Ay(2)-Ay(1)).*u-center(2)).^2).*sqrt((Ax(2)-Ax(1)).^2+(Ay(2)-Ay(1)).^2),0,1);

            end
            
        end
        %%==========================================================
        %%%%%%%%%================================================================
%         
                figure(fig);
                set(fig1,'visible','off');
                set(fig2,'visible','off');
                drawnow;
                fig1=plot(Ax,Ay,'r.-');%,Ax(1),Ay(1),'ro',Ax(2),Ay(2),'r*');
                if(theta==0 |theta==pi |theta==pi/2)
                    fig2=plot((index(:,2)-0.5)*dz,(index(:,1)-0.5)*dz,'bo');
        
                else
                    fig2=plot((index(:,1)-0.5)*dz,(index(:,2)-0.5)*dz,'bo',Q(:,1),Q(:,2),'g.');
                end
                axis(omega);
                  pause;
        %%%%%%%%%================================================================
    end
    if(theta<pi/2)
        R=fliplr(R);
        Rdis=fliplr(Rdis);
    end
    %%%%%%%%%================================================================
    RadonR1=radon(I,theta/pi*180-90);
%     in1=find(R);
%     in2=find(RadonR1);
%     Rref=RadonR1(in2);
%     R1=[Rref(1:2)' R Rref(end-1:end)'];
%     Rspace=linspace(Ttol(in1(1)),Ttol(in1(end)),length(in2));
%     Rspace1=linspace(Ttol(in1(1)),Ttol(in1(end)),length(R1));
    Rspace=linspace(Ttol(1),Ttol(end),length(RadonR1));
    yiC = interp1(Ttol,R,Rspace);
    yiD = interp1(Ttol,Rdis,Rspace);
    errR_C(n)=norm(RadonR1'-yiC);
    errR_D(n)=norm(RadonR1'-yiD);
    errC_D(n)=norm(Rdis-R);
    %%%%%%%%================================================================
    figure('name',sprintf('Projection at angle %d',thetan(n)));

     plot(Ttol,R,'r.-',Ttol,Rdis,'g*-',Rspace,RadonR1,'bo-');
%     plot(Ttol(in1),R(in1),'r.-',Ttol(in1),Rdis(in1),'g*-',Rspace,RadonR1(in2),'bo-');
    title('Red:Continuous; Green: Discrete; Blue: Matlab Radon')
    pause
    %%%%%%%%%================================================================
 
    DisR(:,n)=Rdis';
    ConR(:,n)=R';
    RadonR(:,n)=RadonR1;
    %     iD=iradon(Rdis',theta/pi*180);
    %     iC=iradon(R',theta/pi*180);
    %     iRa=iradon(RadonR(:,n),theta/pi*180);
    %     err(n)=norm(iC-iRa(2:end-1,2:end-1));
    %     figure('name',sprintf('Discrete Reconstruction %d',thetan(n)));imagesc(iD);
    %     figure('name',sprintf('Continuous Reconstruction %d',thetan(n)));imagesc(iC);
    %     figure('name',sprintf('Matlab Reconstruction %d',thetan(n)));imagesc(iRa);
    %     pause;
    
    
end
% Rdis
figure('name','Difference of Radon and MyRadon')
plot(thetan,errR_D,'r.-');%,thetan,errR_C,'b*-',thetan,errC_D,'go-');
% title('Red: Matlab vs. Discrete; Blue: Matlab vs.Continuous; Green: Discrete vs. Continuous')
% figure('name','Discrete Projection');imagesc(DisR);axis xy
% figure('name','Continuous Projection');imagesc(ConR);axis xy
% figure('name','Matlab Radon Projection');imagesc(RadonR);axis xy


time_elapsed=cputime-time;
fprintf('time elapsed: %d \n',time_elapsed);
iDisR=iradon(fliplr(DisR),thetan);
iConR=iradon(fliplr(ConR),thetan);
iRadonR=iradon(RadonR,thetan);
figure('name','Discrete Reconstruction');imagesc(iDisR);axis xy
figure('name','Continuous Reconstruction');imagesc(iConR);axis xy
figure('name','Matlab Reconstruction');imagesc(iRadonR);axis xy

%%% 100x100 with 90 rotations take 236.7500 secs

