N_test=[10:10:2000];
numThetan=1;
n_delta=2*numThetan;
synthetic=1;
for N_ind=1:length(N_test)
    N=N_test(N_ind);
    cross_ind=floor(N/2);
    ex_cor;
    % AnalyticShiftCOR;
    Mt=-log(DisR_Simulated./I0');
    fctn_COR=@(x)sfun_cor(x,full(Mt),sparse(L));% on attenuation coefficients miu;
    tic;
    [f,g]=feval(fctn_COR,[zeros(numThetan*2,1);W(:)]);
    t(N_ind)=toc;
    % err_test{N_ind}=err;
end
save('timeComplexity.mat','N_test','t');
% save('resEffect_1_high.mat','N_test','err_test');
