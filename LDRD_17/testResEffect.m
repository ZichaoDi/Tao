N_test=[2:15]*10;
for N_ind=1:length(N_test)
    N=N_test(N_ind);
    cross_ind=floor(N/2);
    plotPert;
    AnalyticShiftCOR;
    err_test(N_ind,:)=err;
end
save err_test err_test
