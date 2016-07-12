% ===============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: PIR,  Parametric Image Registration
%
%   - data                 Hand, Omega=(0,20)x(0,25), level=5, m=[32,32]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     rotation2D
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all,
viewImage('reset','viewImage','viewImage2D','colormap',gray(256));%,'axis','off'
global testpad DrawError Ro
% ===============================================================================
% T=double(imread('LenaRotate_P30.tiff'));
%T=double(imread('LenaCroppedRotate.tiff'));
% T=double(imread('LenaEye.tiff'));
% R=double(imread('LenaReference.tiff'));
randn('state',0);
for i=72;%:2:89
R=Tiff(['img_2xfm_01',num2str(i),'.h5_H_dpc_cfg.tif'],'r');
T=Tiff(['img_2xfm_01',num2str(i+1),'.h5_H_dpc_cfg.tif'],'r');
T=T.read();
R=R.read();
T=T(1:5:46,1:100:end);
R=R(1:5:46,1:100:end);
figure, subplot(1,2,1),imagesc(R)
subplot(1,2,2),imagesc(T)
pause;
end

c=[3 7]*1;
rec=[10 10]*1;
r=3;
% [T,R]=CreateSimpleTests(c,rec,r);
To=T;
Ro=R;
% Crop=[100 100 50 50];
% T=double(imcrop(R,Crop));

inter('reset','inter','linearInter','regularizer','moments','theta',1e-1);
omega=[0 size(R,1) 0 size(R,2)];
omegat=omega;%[10 10+size(T,1) 40 40+size(T,2)];%[Crop(1)-1 Crop(1)+Crop(3) Crop(2)-1 Crop(2)+Crop(4)];%
mt=floor(size(T)/1);
m = floor(size(R)/1);
T = inter('coefficients',T,[],omegat,'regularizer','moments','theta',1e-1);
R = inter('coefficients',R,[],omega,'regularizer','moments','theta',1e-1);
%vT=10*randn(size(T));
%T=T+vT;
%vR=10*randn(size(R));
% R=R+vR;
% Rtt=R(omegat(1):omegat(2)-1,omegat(3):omegat(4)-1);
% ===============================================================================
testpad=0;%mean(mean(R));
DrawError=0;
distance('reset','distance','SSD');
center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','translation2D','c',center);
w0 = trafo('w0');
beta = 0; M =[]; wRef = []; % disable regularization
% beta = 1; M = 5e2*speye(length(w0)); wRef = w0;     % initilize the regularization

% initialize plots
FAIRplots('reset','mode','PIR-GN','fig',1,'scale',0);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'omegat',omegat,'m',m,'mt',mt));
pause;
xc = getCellCenteredGrid(omega,m);
Rc = inter(R,omega,xc);
% figure(112);viewImage(Rc,omega,m,'scale',1);
% %hold on; plotGrid(xc,omega,m,'spacing',ceil(m/32));axis(omega);hold off;
% hold on; viewImage(T(:),omegat,mt,'scale',1);
% pause;
xt=getCellCenteredGrid(omegat,mt);
% fctn = @(wc) PIRobjFctn1(T,Rc,omega,omegat,m,mt,beta,M,wRef,xc,xt,wc);
% fctn = @(wc) PIRobjFctnR(T,Rc,omega,omegat,m,mt,beta,M,wRef,xc,xt,wc);
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc);
%==========================================================================
% disP=c(1)-c(2);
% w=linspace(-5,1,rec(1)*10);
% Jc=[];
% for i=1:length(w)
%     Jc(i)=feval(fctn,[w(i),w(i)]');
%     
% end
% drawnow;
% figure(12);plot(w,Jc,'r.-');hold on;
% return;
% % optimize
% %  checkDerivative(fctn,w0);
% %  foo(fctn,w0);
% %  return;


[wc,his] = GaussNewton(fctn,w0,'Plots',@FAIRplots,'solver','direct','yStop',w0);
% [wc,his] = lBFGS(fctn,w0,'Plots',@FAIRplots,'solver','direct');



