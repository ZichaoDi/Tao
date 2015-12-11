N=20;
numThetan=25;
xstar=exAlternate(N,numThetan);
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
