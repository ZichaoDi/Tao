%%% Simulate XRF of a given fixed object with rotating predifined detector and beam
%%% Travelling of fluorescence photon is approximated as the area covered by solid angle 
global BeforeEmit plotTravel SSDlet area_xrf
% global EmptyBeam RealBeam
global LogScale Tol
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
L=sparse(zeros(numThetan*(nTau+1),mtol));%zeros(numThetan*(nTau+1)*prod(m),1);
GlobalInd=cell(numThetan,nTau+1);
area_xrf=sparse(zeros(numThetan*(nTau+1),mtol));
SelfInd=repmat({cell(3,1)},[numThetan*(nTau+1),mtol]);

fprintf(1,'====== Fluorescence Detector Resolution is %d\n',numChannel);

for n=1:numThetan
    theta=thetan(n)/180*pi;
    fprintf(1,'====== Angle Number  %d of %d: %d\n',n,numThetan,thetan(n));
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    DetKnot=DetKnot0*TransMatrix;
    SourceKnot=SourceKnot0*TransMatrix;
    SSDknot=SSDlet*TransMatrix;
    %=================================================================
    for i=1:nTau+1 %%%%%%%%%================================================================
        % Initialize
        xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
        ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
        BeforeEmit=1;
        %=================================================================
        [index,Lvec,linearInd]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta);
        ind_bt=(i-1)*numThetan+n;
        GlobalInd{ind_bt}=index;
        for j=1:size(index,1)
            CurrentCellCenter=[(index(j,1)-1/2)*dz(1)-abs(omega(1)),(index(j,2)-1/2)*dz(2)-abs(omega(3))];
            currentInd=sub2ind(m,index(j,2),index(j,1));
            L(ind_bt,currentInd)=Lvec(j);
            if(j>1 && j<size(index,1))
                SelfInd{ind_bt,currentInd}{1}=sub2ind(m,index(j-1:-1:1,2),index(j-1:-1:1,1));
                SelfInd{ind_bt,currentInd}{3}=kron(Lvec(j-1:-1:1),MU_e(:,1,1)');
            else
                SelfInd{ind_bt,currentInd}{1}=sub2ind(m,index(j-1:-1:1,2),index(j-1:-1:1,1));
                SelfInd{ind_bt,currentInd}{3}=kron(Lvec(j-1:-1:1),MU_e(:,1,1)');
            end
            %% ===========================================================
            %% Self-absorption
            in_after=find(inpolygon(xc(:,1),xc(:,2),[CurrentCellCenter(1) SSDknot(1,1) SSDknot(NumSSDlet,1)],[CurrentCellCenter(2) SSDknot(1,2) SSDknot(NumSSDlet,2)])); 
            in_after=setdiff(in_after,currentInd); %% energy is not attenuated from the source point
            area_xrf(ind_bt,currentInd)=prod(dz)*length(in_after);
            SelfInd{ind_bt,currentInd}{2}=in_after;
        end
    end
end

clear Wsub A knot DetKnot DetKnot0 SSDknot SSDlet SourceKnot SourceKnot0 MU_after
clear xbox ybox eX eY  x y xc 




