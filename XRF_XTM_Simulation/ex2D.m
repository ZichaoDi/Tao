%%%Simulate XRF of a given object with predifined detector and beam
do_setup;
optXRF;
Joint=1;
optXRF;
Alternate=0;
optXRF;
Joint=0;
optXRF;
save('drawing.mat','icycle','err0','errOut_XRF','errOut_Joint','maxOut','ErrIter_XRF','ErrIter_Joint');
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
