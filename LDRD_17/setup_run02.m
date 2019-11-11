nslice=1200;
load(['Run02_',num2str(slice)]);
XRF_decom=-log(data(1:tau_rate:end,1:ang_rate:end)'./max(max(data)));
nTau=size(XRF_decom,2)-1;
N=nTau+1;
numThetan=size(XRF_decom,1);
thetan_real=linspace(0,180,numThetan);
Z=20;
