
clear, close all,
%help(mfilename);
viewImage('reset','viewImage','viewImage2D','colormap',gray(256));%,'axis','off'
global testpad xc TbRs
randn('state',0);

% ===============================================================================

T=double(imread('LoganReference.tiff'));
R=flipud(double(imread('LoganTest1.tiff')));
%  R=double(imread('CroppedRotate.tiff'));
TbRs=1;
% Crop=[100 100 50 50];
% R=double(imcrop(T,Crop));
% omega=[Crop(1)-1 Crop(1)+Crop(3) Crop(2)-1 Crop(2)+Crop(4)];
% omega=10+omega;

%%% ===========================================
% imshow(uint8(R))
% figure,imshow(uint8(T))
%   return;
inter('reset','inter','splineInter','regularizer','moments','theta',1e-1);
Inix=20;
Iniy=30;
omega=[Inix Inix+size(R,1) Iniy Iniy+size(R,2)];
omegat=[0 0+size(T,1) 0 0+size(T,2)];%omega;%
mt=floor(size(T)/1);

m = floor(size(R)/1);
T = inter('coefficients',T,[],omegat,'regularizer','moments','theta',1e0);
R = inter('coefficients',R,[],omega,'regularizer','moments','theta',1e0);
% vT=1*rand(size(T));
% T=T+vT;
% vR=10*rand(size(R));
% R=R+vR;
% ===============================================================================
testpad=0;%mean(mean(R));
distance('reset','distance','SSD');

center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rigid2D','c',center);
w0 = trafo('w0');
beta = 0; M =[]; wRef = []; % disable regularization
% beta = 1e0; M = 5e2*speye(length(w0)); wRef =[0.4,0,0]';     % initilize the regularization
hd=prod(1./m);
% beta = 1e0; M = 5e0*hd*spdiags(ones(2,1)*[-1,2,1],-1:1,2,2); wRef = w0;%[-5,-5]';

% beta = 1e0; M = 1e2*hd*spdiags(ones(2,1)*[1,-1],0:1,2,2); wRef =w0;%[-15,-15]';

%====== initialize plots
FAIRplots2('reset','mode','PIR-GN','fig',1,'scale',1);
FAIRplots2('init',struct('Tc',T,'Rc',R,'omega',omega,'omegat',omegat,'m',m,'mt',mt));
pause;
xc = getCellCenteredGrid(omega,m);
Rc = inter(R,omega,xc);
xt=getCellCenteredGrid(omegat,mt);
% figure(11);
% surf(reshape(xc(1:length(xc)/2),m(1),m(2)),reshape(xc(length(xc)/2+1:end),m(1),m(2)),R);
% figure,
% surf(reshape(xt(1:length(xt)/2),mt(1),mt(2)),reshape(xt(length(xt)/2+1:end),mt(1),mt(2)),T);
% return;
fctn = @(wc) PIRobjFctnTbRs(T,Rc,omega,omegat,m,mt,beta,M,wRef,xc,xt,wc);% local registration
% fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc); %global registration


% w2=linspace(0,50,10);
% w3=linspace(0,20,10);
% [W2,W3]=meshgrid(w2,w3);
% W2=W2(:);
% W3=W3(:);
% W=sqrt(W2.^2+W3.^2);
% n=length(W);
% w1=linspace(-1,1,prod(size(W2)));
% Jc=zeros(size(w1));
% for i=1:length(W2)
%     Jc(i)=feval(fctn,[0.52,W2(i),W3(i)]');
% end
% figure(12);
% surf(reshape(W2,sqrt(n),sqrt(n)),reshape(W3,sqrt(n),sqrt(n)),reshape(Jc,sqrt(n),sqrt(n)));
% return;
% w=linspace(-15,10,100);
% Jc=zeros(size(w));
% for i=1:length(w)
%     Jc(i)=feval(fctn,[0.52,w(i),w(i)]');
% end
% drawnow;
% figure(12);plot(w,Jc,'r.-');hold on;
% return;
% foo(fctn,w0);
% checkDerivative(fctn,w0)
% return;

%===== optimize
   w0=[0.25 0 0]';
[wc,his] = GaussNewton(fctn,w0,'Plots',@FAIRplots2,'solver','direct','yStop',w0);
% [wc,his] = lBFGS(fctn,w0,'Plots',@FAIRplots2,'solver','direct','yStop',w0);




