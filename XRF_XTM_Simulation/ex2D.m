do_setup;
% Joint=1;
%optXTM_XRF;
% Joint=0;
% optXTM_XRF;
% return;
%mg;
%for i=1:10;
%    mgit;
%end

opt;
init=1;
ML;
init=0;
for nn=1:25; ML;end
    
save(['opt_ML',num2str(N(1)),'_',num2str(nTau+1),'_',num2str(numThetan),'GS_',sample,'.mat'],'per_rec','sing');

