% function xstar=exAlternate(simulate)
%%%Simulate XRF of a given object with predifined detector and beam
do_setup;
shift_test=[-1:5:40];
x_shift=zeros(length(shift_test),N^2*NumElement);
for sind=1:length(shift_test)
    load spectra_30_aligned;
    spectra=0.*spectra_30_aligned;
    center_shift=shift_test(sind);
    if(center_shift>=0)
        spectra(center_shift:end,:,:)=spectra_30_aligned(1:end-center_shift+1,:,:);
        spectra(1:center_shift-1,:,:)=spectra_30_aligned(end-center_shift+2:end,:,:);
    else
        spectra(1:end+center_shift,:,:)=spectra_30_aligned(abs(center_shift)+1:end,:,:);
        spectra(end+center_shift+1:end,:,:)=spectra_30_aligned(1:abs(center_shift),:,:);
    end

    data_xrf_raw=permute(spectra(1:tau_rate:end,1:ang_rate:end,:),[3 2 1]);
    data_xrf_raw=sparse(reshape(double(data_xrf_raw),[size(data_xrf_raw,1),size(data_xrf_raw,2)*size(data_xrf_raw,3)]));
    XRF_raw=data_xrf_raw';
    DecomposedElement=0;
    opt;
    x_shift(sind,:)=xstar';
    save x_shift x_shift
end
% if(simulate)
%     do_setup_temp;
%     DecomposedElement=0;
%     opt_temp;
% else
%     do_setup;
%     DecomposedElement=0;
%     Be=[0 1 2 3 4 8 16];
%     TempBeta=1;
%     for j_iter=1:length(Be);
%         Beta=Be(j_iter);
%         opt;
%     end
%     TempBeta=0; Beta=1;
%     opt;
%     % Joint = 0;
%     % opt;
%     % maxiter =100;
%     % Joint = 1;
%     % opt;
% end
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
