clear, close all, help(mfilename);

% load data, set viewer, interpolator, transformation, distance
%dataT=imread('LenaCroppedRotate.tiff');
dataT=imread('lenaEye.tiff');
dataR=imread('LenaReference.tiff');
omega=[0 size(dataR,1) 0 size(dataR,2)];
m = floor(size(dataR)/4);
viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
inter('reset','inter','linearInter');
distance('reset','distance','SSD');
center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rigid2D','c',center);

% discretize T and R, make shortcut do difference image
xc = getCellCenteredGrid(omega,m);
Tc = inter(dataT,omega,xc);
Rc = inter(dataR,omega,xc);
Dc = @(Tc) 128+(Tc-Rc)/2;
hd = prod((omega(2:2:end)-omega(1:2:end))./m);

% choices of coefficients
tt=[0     0     0;
    -0.2  0     0;
    -0.2  -10   10;
    0     0     0;
    0.1   0     0;
    0.1   10    10;
    0.1   100   -10;
    0.1   150   -20;
    0.1   200   -30;
    0.1   260   -7; 
    0.4   0     0;
    0.4   10    10;
    0.4   100   -10;
    0.4   150   -20
    0.4   200   -30; 
    0.4   260   -7 
    0.467 260   -8
    0.47  260   -8.5
    0.524 260   -8.5
    0.524 200   10;
    0.524 -25   -250
    0.524 -10   -250
    0.524 0     -250
    0.524 40   -200]';
ntot=size(tt,2);
D  = zeros(ntot,1);
for k=1:ntot,
    yc   = trafo(tt(:,k),xc);         % transform the grid X -> Y
    Tc   = inter(dataT,omega,yc);     % compute transformed image
    D(k) = hd*norm(Tc-Rc)^2;
    % visualize results
    if k == 1,
        FAIRfigure(1,'figname',mfilename);
        subplot(2,2,1); viewImage(Rc,omega,m);
        title('reference','fontsize',30)
        subplot(2,2,2); th = viewImage(Tc,omega,m);
        title('template','fontsize',30)
        subplot(2,2,3); dh = viewImage(Dc(Tc),omega,m);
        title('difference','fontsize',30)
        subplot(2,2,4);
        hold on;
        ph= plot(1:k,D(1:k),'w.-','markersize',20);
        title('SSD versus transformation','fontsize',10)
        drawnow;
        
    else
        set(dh,'cdata',reshape(Dc(Tc),m)')
        subplot(2,2,4);
        
        set(ph,'visible','off');
        ph=  plot(1:k,D(1:k),'k.-','markersize',20);
        drawnow;
    end;
    pause;
end;
