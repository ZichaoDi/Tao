do_setup;
 Joint=1;
optXTM_XRF;
 Joint=0;
 optXTM_XRF;
 return;
mg;
for i=1:10;
    mgit;
end
save(['xs_mg',num2str(N(1)),'_',num2str(nTau+1),'_',num2str(numThetan),'TNBC_',sample,'.mat'],'v0','v');
% profile on;
% p = profile('info');
% save myprofiledata9 p
% clear p
% load myprofiledata
% profview(0,p)


