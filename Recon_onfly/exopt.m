beta=[1 0.0; 1 0.5; 1 5; 1 50; 1 100; 1 1e3];
for t=1:size(beta,1)
    TempBeta=beta(t,1); Beta=beta(t,2);
    if(t==1) 
        initialize=1; 
    else 
        initialize=0;
    end
    opt;
end
