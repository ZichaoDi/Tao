function xstar=exAlternate(simulate)
%%%Simulate XRF of a given object with predifined detector and beam
if(simulate)
    do_setup_temp;
    DecomposedElement=0;
    opt_temp;
else
    do_setup;
    DecomposedElement=0;
    Beta=1;
    for bt=1:5;
        opt;
        Beta = Beta/2;
        Beta=5;
    end
    % Joint = 0;
    % opt;
    % maxiter =100;
    % Joint = 1;
    % opt;
end
%{
return;
save(['drawing2',num2str(N(1)),'_',num2str(numThetan),'.mat'],'errOut_Joint');
Joint=0;
optXRF;
save(['drawing1',num2str(N(1)),'_',num2str(numThetan),'.mat'],'icycle','maxOut','err0','errOut_XRF');
Alternate=0;
Joint=1;
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
