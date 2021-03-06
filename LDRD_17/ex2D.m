%%%Simulate XRF of a given object with predifined detector and beam
x_res=[];
f_dr=[];
N_total=[13 25 50 60 70 80 90 100 120 140];% 17 9];%[129 65  9 5];%
res=2;
for i_N=1:length(N_total);
    N=N_total(i_N);
    plotPert;
    opt;
    f_dr(i_N,1)=f;
    x0_opt=x0;
    opt_dr;
    f_dr(i_N,2)=f;
end
% save(['x_res_cP_positive',sample,'.mat'],'x_res');
% exAlternate(0);
%{
return;
h_admm=[];for i=1:size(x_admm,2)-1, h_admm(:,i)=x_admm(:,i+1)-x_admm(:,i);end
f1=[];for i=1:size(x_admm,2),f1(i)=feval(fctn,x_admm(:,i));end
gg=[];for i=1:size(x_admm,2)-1;[~,gg(:,i)]=feval(fctn1,x_admm(:,i));end
 res=f1(2:end)-f1(1:end-1);
 for i=1:size(x_admm,2)-1;df(i)=h_admm(:,i)'*gg(:,i);end
     save('linearizePer.m');
     figure, semilogy(1:size(x_admm,2)-1,abs(df),'r.-',1:size(x_admm,2)-1,abs(res),'b.-');
     return;
save(['drawing2',num2str(N(1)),'_',num2str(numThetan),'.mat'],'errOut_Joint');
Joint=0;
optXRF;
save(['drawing1',num2str(N(1)),'_',num2str(numThetan),'.mat'],'icycle','maxOut','err0','errOut_XRF');
Alternate=0;
optXRF;
save(['drawing3',num2str(N(1)),'_',num2str(numThetan),'.mat'],'ErrIter_Joint');
Joint=0;
optXRF;
save(['drawing4',num2str(N(1)),'_',num2str(numThetan),'.mat'],'ErrIter_XRF');
return;
startup;
close all;
global plotTravel plotSpec plotWhole
plotTravel=0; % If plot the intersection of beam with object
plotSpec = 1; % Do you want to see the spectra? If so plotSpec = 1
plotWhole =0;
if plotSpec
    fig1=[]; fig2=[]; fig3=[];
end
DefineObject; %% Produce W, MU
Define_Detector_Beam; %% provide the beam source and Detectorlet
UnitSpectrum; %% Produce BindingE M
theta=1:10:90;% Projection Angles
XRF=SimulateXRF(W,MU,BindingE,M,theta,DetChannel, numChannel, nTau, DetKnot0, SourceKnot0);
%}
