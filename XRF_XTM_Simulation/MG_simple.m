%%% simple multigrid
global x0 current_n N;
global Beta;

close all;
more on;
plotResult=0;
do_setup;
tic;
cycle=1;
maxCycle=1;
fid = fopen('betaPareto_Tensor.txt','a');
    
while(cycle<=maxCycle)

for it=-6:4
Beta=10^(it);    

for level=1:length(N)-1
    current_n=N(level);
    maxiter=1;
    if(cycle>1 & current_n==N(1))
    x0=xstar;
    end
    fprintf('===================== %d\n',current_n)
    optXTM_XRF_Tensor;
    if(current_n~=N(end))
        x0=downdate(xstar,1);
    end
end
for level=length(N):-1:1
    current_n=N(level);
    if(current_n==N(end))
        if(length(N)==1)
            maxiter=10000;
        else
        maxiter=100;
        end
    else
        disp('====================== Start post-smoothing')
        maxiter=1;
        if(current_n==N(1))
            maxiter=100;
        end
    end
    if(current_n==N(1))
        plotResult=1;
    else
        plotResult=0;
    end
    fprintf('===================== %d\n',current_n)
    optXTM_XRF_Tensor;
    if(current_n~=N(1))
        x0=update(xstar,1);
    end
end
cycle=cycle+1;
err_h=norm(xstar-W(:));
T=toc;
report_results(N);

[f,g,f1,f2]=feval(fctn,xstar);
fprintf(fid,'%12.4e    %12.4e     %12.4e    %12.4e      %12.4e     %f\n',Beta,f,f1,f2,err_h,T);
end

end
  fclose(fid);
