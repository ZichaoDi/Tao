% beta=[1 1; 1 10; 1 20; 1 30;1 40];%[1 0; 0 1;1 1; 1 10; 1 15; 1 20; 1 30] ;%0.5;1 1; 1 10; 1 100; 1 1e3];
initialize=1; 
for slice =[30 40 50 60];
    beta=[1 0; 1 1];% 1 30; 1 40; 1 50; 1 60; 1 70; 1 80];
    for t=1:size(beta,1)
        TempBeta=beta(t,1); Beta=beta(t,2);
        opt;
    end
    initialize=0;
end
