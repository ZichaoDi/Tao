beta=[1 0.05; 1 0.01; 1 0.005;; 1 0.001];
for t=1:size(beta,1)
    TempBeta=beta(t,1); Beta=beta(t,2);
    if(t==1) 
        initialize=1; 
    else 
        initialize=0;
    end
    opt;
end
