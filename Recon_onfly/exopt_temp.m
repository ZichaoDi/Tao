% beta=[1 0.0; 1 0.05; 1 0.1; 1 0.25; 1 0.5; 1 1; 1 2; 1 4];
beta=[1 0; 1 0.25; 1 0.5;  1 1; 1 2; 1 4];
for t=1:size(beta,1)
    TempBeta=beta(t,1); Beta=beta(t,2);
    if(t==1) 
        initialize=1; 
    else 
        initialize=0;
    end
    opt_temp;
end
