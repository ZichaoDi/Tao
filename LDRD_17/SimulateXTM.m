%%% Simulate X-ray Transmission Microscopy
close all
%%%%%%%===================== Set up experiments
startup;
plotProjection=0;
plotTravel=0; % If plot the intersection of beam with object
plotSpec = 0; % Do you want to see the spectra? If so plotSpec = 1
plotUnit=0;
plotSpecSingle=0;
NoSelfAbsorption=0;
PlotObject=1;
more off;

Define_Detector_Beam_Gaussian; %% provide the beam source and Detectorlet
UnitSpectrumSherman_Gaussian; %% Produce BindingEnergy M
DefineObject_Gaussian; %% Produce W, MU
thetan=[1:15:180];% Projection Angles
%%%%%%%==============================================================
%%%=========Start Simulation

time=cputime;
e=ones(m(1),1);
DisR=[];
RadonR=[];
Ltol=cell(length(thetan),nTau);
fig1=[]; fig2=[]; fig3=[];
for n=1:length(thetan)
    
    fprintf(1,'====== Angle Number  %d of %d\n',n,length(thetan));
    
    theta=thetan(n)/180*pi;
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    DetKnot=DetKnot0*TransMatrix;
    SourceKnot=SourceKnot0*TransMatrix;
    if(plotTravel)
        figure(fig);
        set(fig3,'visible','off');
        for i=1:nTau+1
            fig3=plot([SourceKnot(i,1),DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)],'r.-');hold on;
        end
        pause;
    end
    for i=1:nTau+1
        %%%%%%%%%================================================================
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        [Ax, Ay] = polyxpoly([double(vpa(SourceKnot(i,1),5)),DetKnot(i,1)],[SourceKnot(i,2),DetKnot(i,2)], xbox, ybox);
        if(isempty(Ax) | (length(Ax)==1 &length(Ay)==1 ))
            %             fprintf('no intersection \n')
            index=[];
            Q=[];
            L=sparse(zeros(m(1),m(2)));
            Rdis(i)=0;
            Ltol{n,i}=[];
        else
            if(theta<=pi/2)
                [Ax,in]=sort(Ax);A=unique([Ax,Ay(in)],'rows');Ax=A(:,1);Ay=A(:,2);
            else
                [Ay,in]=sort(Ay);[A,in]=unique([Ax(in),Ay],'rows');
                if(size(A,1)>1)
                    A=A(in,:);Ax=A(:,1);Ay=A(:,2);
                end
            end
            %%%%%%%%%================================================================
            if(theta==pi/2)
                Q=[repmat(Ax(1),size(y')),y'];
            elseif(theta==0 | theta==pi)
                Q=[x',repmat(Ay(1),size(x'))];
            else
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
                index=[floor(myvpa((Q(:,1)+abs(omega(1)))/dz))+1,floor(myvpa((Q(:,2)+abs(omega(3)))/dz))+1];
                
            else
                index=[floor(myvpa((Q(:,1)+abs(omega(1)))/dz))+1,ceil(myvpa((Q(:,2)+abs(omega(3)))/dz))];
                
            end
            index=index(find(index(:,1)>0 & index(:,1)<=m(1)& index(:,2)<=m(1) & index(:,2)>0),:);
            %%%%%%%%================================================================
            if(plotTravel)
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
                end
            end
            %%%%%%%%%================================================================
            L=sparse(zeros(m(1),m(2)));
            for j=1:size(index,1)
                L(index(j,2),index(j,1))=Lvec(j);
            end
            
            Ltol{n,i}=L;
            Rdis(i)=I0*exp(-e'*(MU.*L)*e); %% Discrete case
        end
        
        
    end
    DisR(:,n)=Rdis';
    %%%%%%%%================================================================
    if(plotProjection)
        Ttol=1:nTau+1;
        figure('name',sprintf('Projection at angle %d',thetan(n)));
        plot(Ttol,Rdis,'g*-');
        pause
    end
end
time_elapsed=cputime-time;


