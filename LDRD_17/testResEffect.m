N_test=[20:10:500];
for N_ind=1:length(N_test)
    N=N_test(N_ind);
    cross_ind=floor(N/2);
    ex_cor;
    AnalyticShiftCOR;
    err_test(N_ind,:)=err;
end
save('err_test500_1.mat','N_test','err_test');
