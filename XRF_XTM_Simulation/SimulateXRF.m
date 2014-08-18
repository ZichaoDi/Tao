%%%Simulate XRF of a given object with predifined detector and beam
%% Disregard self-absorption
% function XRF=SimulateXRF(W,MU,BindingEenergy,M,thetan,DetChannel, numChannel, nTau, DetKnot0, SourceKnot0);
global plotSpecSingle
close all;
plotTravel=1; % If plot the intersection of beam with object
plotSpec = 0; % Do you want to see the spectra? If so plotSpec = 1
plotWhole =0;
plotUnit=0;
plotSpecSingle=0;
if plotTravel
    fig1=[]; fig2=[]; fig3=[];
end
% global m x y omega NumElement plotTravel plotSpec plotWhole dz
startup;
more off;
DefineObject; %% Produce W, MU
Define_Detector_Beam; %% provide the beam source and Detectorlet
% UnitSpectrum; %% Produce BindingEnergy M
UnitSpectrumSherman
thetan=[30 120 210 300];% Projection Angles
%%%%%%%==============================================================
%%%=========Start Simulation
Energy=BindingEnergy;

XRF=cell(length(thetan),nTau+1);

Ltol=cell(length(thetan),nTau+1);
GlobalInd=cell(length(thetan),nTau+1);
RMlocal=zeros(m(1),m(2),NumElement);
fprintf(1,'====== Detector Resolution is %d\n',nTau);
for n=1:length(thetan)
    
    fprintf(1,'====== Angle Number  %d of %d\n',n,length(thetan));
    
    theta=thetan(n)/180*pi;
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    DetKnot=DetKnot0*TransMatrix;
    SourceKnot=SourceKnot0*TransMatrix;
    if plotWhole
        figure(fig);
        set(fig3,'visible','off');
        for i=1:nTau+1
            fig3=plot([SourceKnot(i,1),DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)],'r.-',SourceKnot(i,1),SourceKnot(i,2),'bo');hold on;
        end
    end
    
    
    for i=1:nTau+1 %%%%%%%%%================================================================
        % Initialize
        XRF{n,i} = zeros(numChannel,1);
        
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        [Ax, Ay] = polyxpoly([SourceKnot(i,1),DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)], xbox, ybox);
        
        if(isempty(Ax) | (length(Ax)==1 &length(Ay)==1 ))
            %             fprintf('no intersection \n')
            index=[];
            Q=[];
            L=sparse(zeros(m(1),m(2)));
            
        else
            if(theta<=pi/2)
                [Ax,in]=sort(Ax);A=unique([Ax,Ay(in)],'rows');Ax=A(:,1);Ay=A(:,2);
            else
                [Ay,in]=sort(Ay);[A,in]=unique([Ax(in),Ay],'rows');
                if(size(A,1)>1)
                    A=A(in,:);Ax=A(:,1);Ay=A(:,2);
                end
            end
            
            if(theta==pi/2)
                Q=[repmat(Ax(1),size(y')),y'];
            elseif(theta==0 | theta==2*pi)
                Q=[x',repmat(Ay(1),size(x'))];
            elseif(theta==pi)
                Q=[x(end:-1:1)',repmat(Ay(1),size(x'))];
                
            elseif(theta==3*pi/2)
                Q=[repmat(Ax(1),size(y')),y(end:-1:1)'];
                
            else
                Q=[[x', (Ay(2)-Ay(1))/(Ax(2)-Ax(1)).*x'+(Ay(1)*Ax(2)-Ay(2)*Ax(1))/(Ax(2)-Ax(1))];...
                    [(y'-(Ay(1)*Ax(2)-Ay(2)*Ax(1))/(Ax(2)-Ax(1)))./((Ay(2)-Ay(1))/(Ax(2)-Ax(1))),y']];
            end
            
            indx=find(Q(:,1)-omega(1)<-1e-6 |Q(:,1)-omega(2)>1e-6); indy=find(Q(:,2)-omega(3)<-1e-6 |Q(:,2)-omega(4)>1e-6);
            Q=setdiff(Q,Q([indx;indy],:),'rows');
            Q=unique(Q,'rows');
            dis=distance(Q,repmat(SourceKnot(i,:),size(Q,1),1));
            [dis,InterOrder]=sort(dis);
            Q=Q(InterOrder,:);
            Lvec=sqrt((Q(2:end,1)-Q(1:end-1,1)).^2+(Q(2:end,2)-Q(1:end-1,2)).^2);
            %%%%%%%%%================================================================
            if(theta>=0 & theta<=pi/2 | theta>=pi & theta<=3*pi/2 )
                index=[floor(myvpa((Q(:,1)+abs(omega(1)))/dz))+1,floor(myvpa((Q(:,2)+abs(omega(3)))/dz))+1];
            else
                index=[floor(myvpa((Q(:,1)+abs(omega(1)))/dz))+1,ceil(myvpa((Q(:,2)+abs(omega(3)))/dz))];
                
            end
            index=index(index(:,1)>0 & index(:,1)<=m(1)& index(:,2)<=m(1) & index(:,2)>0,:);
            [temp,subInd]=unique(index,'rows');
            clear temp;
            index=index(sort(subInd),:);
            GlobalInd{n,i}=index;
            %%%%%%%%%================================================================
            if plotTravel
                figure(fig);
                set(fig1,'visible','off');
                set(fig2,'visible','off');
                drawnow;
                fig1=plot(Ax,Ay,'r.-',SourceKnot(i,1),SourceKnot(i,2),'ro',DetKnot(i,1),DetKnot(i,2),'r*');
                if(~isempty(index))
                    %                     fig2=plot((index(:,1)-abs(omega(1)))*dz,(index(:,2)-abs(omega(3)))*dz,'bo',Q(:,1),Q(:,2),'g.');
                    fig2=plot((index(:,1)-1/2)*dz-abs(omega(1)),(index(:,2)-1/2)*dz-abs(omega(3)),'bo',Q(:,1),Q(:,2),'g.');
                    
                end
                pause;
            end
            %%%%%%%%%================================================================
            
            L=sparse(zeros(m(1),m(2)));
            RM=cell(m(1),m(2));
            xrfSub=zeros(size(BindingEnergy));
            for j=1:size(index,1)
                L(index(j,2),index(j,1))=Lvec(j);
                if(j==1)
                    I_incident=1;
                    temp_sum=0;
                else
                    temp_sum=temp_sum+L(index(j-1,2),index(j-1,1))*MU(index(j-1,2),index(j-1,1));
                    I_incident=exp(-temp_sum);
                end
                %% Self-absorption
                %                 temp_after=0;
                %                 for t_after=index(j,1)+1:m(2)
                %                     temp_after=temp_after+dz*MU_after(index(j,2),t_after);
                %                 end
                %                 I_after=exp(-temp_after);
                %% ===========================================================
                Wsub=reshape(W(index(j,2),index(j,1),:),size(M));
                RM{index(j,2),index(j,1)}=Lvec(j)*I_incident*(Wsub.*M);%*I_after;
                temp=zeros(size(RMlocal(index(j,2),index(j,1),:)));
                temp(1,1,:)=Lvec(j)*I_incident*(Wsub.*M);
                RMlocal(index(j,2),index(j,1),:)=RMlocal(index(j,2),index(j,1),:)+temp;
                ExistInd=find(Wsub~=0);
                for tsub=1:length(ExistInd)
                    DetInd=find(Energy(ExistInd(tsub))==BindingEnergy);
                    
                    xrfSub(DetInd)=xrfSub(DetInd)+RM{index(j,2),index(j,1)}(ExistInd(tsub));
                end
            end
            Ltol{n,i}=L;
            p1=GaussianFit1(BindingEnergy,xrfSub);
            if(plotSpec)
                finalfig=figure('name',sprintf('XRF at beam %d with angle %d',i,thetan(n)));
                plot(DetChannel,p1,'r-')
                xlabel('Energy Channel','fontsize',12); ylabel('Intensity','fontsize',12)
                pause;
            end
            XRF{n,i}=p1;
        end
        
    end
end

