global reg lambda n_delta initialize
global f1 f2
ReconAttenu=1;
N=128;
numThetan=40;
synthetic=1;
n_delta=2*numThetan;
samSet={'MRI','Phantom','Golosio'};
prj=zeros(numThetan,N,length(samSet));
wtrue=zeros(N,N,length(samSet));
for isample=1:length(samSet)
    sample=samSet{isample};
    Ntest=2;
    DisR=[];
    L_cr=[];
    initialize=1;
    ind_scan=1;
    for ind_cr=1:Ntest;
            do_setup; 
            L_cr(:,:,ind_cr)=L;
            DisR(:,:,ind_cr)=DisR_Simulated;
    end
    return;
    prj(:,:,isample)=-log(DisR(:,:,2)./I0);
    wtrue(:,:,isample)=W;
end
% save('syntheticMultiSino.mat','prj','thetan','wtrue');
return;
% save('data30.dat','b','-ASCII');
lam=-2;%[-15:0, 0.2 0.4 0.8 1.6];
reg_str={'TV','L2'};
for regi=1
    if(regi==1)
        maxiter=10;
    else
        maxiter=100;
    end
    reg=reg_str{regi};
for i=1:length(lam)
    lambda=0*10^lam(i);
    opt;
    xreg(:,i)=x;
    ftotal(i)=f;
    fmis(i)=f1;
    freg(i)=f2;
    greg(i)=norm(g);
end
% save(['Phantom_1cor_init_',reg,num2str(lam(end)),'.mat'],'xreg','ftotal','fmis','freg','greg');
end
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
