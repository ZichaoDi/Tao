global initialize
numThetan=1;
n_delta=2*1;
% if(n_delta==2)
    Ntest=2;
    DisR=[];
    L_cr=[];
    nslice=150;
    initialize=1;
    ind_scan=1;
    for ind_cr=1:Ntest;
            do_setup; 
            L_cr(:,:,ind_cr)=L;
            DisR(:,:,ind_cr)=DisR_Simulated;
    end

%     Mt1=-log(DisR(:,:,1)./I0);
%     Mt2=-log(DisR(:,:,2)./I0);
%     return;
%     opt;
%     % save([sample,'_1cor_ini_explicit',num2str(numThetan),'.mat'],'N','W','aligned','deltaStar','f_global','nTau','n_delta','nit_res','numThetan','x0_opt','x_res','x','delta0','dz');
%     opt_dr;
%     % save('iter_init_1cor.mat','f0','f1','f1_dr','xiter0','xiter1','x1_dr');
%     % save([sample,'_1cor_ini_implicit',num2str(numThetan),'.mat'],'N','W','aligned','deltaStar','f_global','nTau','n_delta','nit_res','numThetan','x0_opt','x_res','x','delta0','dz');
%     %%=================================================
%     n_delta=2*numThetan;
% % elseif(n_delta==2*numThetan)
%     Ntest=2;
%     DisR=[];
%     L_cr=[];
%     nslice=150;
%     initialize=1;
%     ind_scan=1;
%     for ind_cr=1:Ntest;
%             do_setup; 
%             L_cr(:,:,ind_cr)=L;
%             DisR(:,:,ind_cr)=DisR_Simulated;
%     end
% 
%     Mt1=-log(DisR(:,:,1)./I0);
%     Mt2=-log(DisR(:,:,2)./I0);
%     opt;
%     % save([sample,'_mcor_ini_explicit',num2str(numThetan),'.mat'],'N','W','aligned','deltaStar','f_global','nTau','n_delta','nit_res','numThetan','x0_opt','x_res','x','delta0','dz');
%     opt_dr;
%     % save('iter_init_mcor.mat','f0','f1','f1_dr','xiter0','xiter1','x1_dr');
%     % save([sample,'_mcor_ini_implicit',num2str(numThetan),'.mat'],'N','W','aligned','deltaStar','f_global','nTau','n_delta','nit_res','numThetan','x0_opt','x_res','x','delta0','dz');
% % end% 
