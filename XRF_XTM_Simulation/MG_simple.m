%%% simple multigrid
global x0 current_n N;
global Beta maxiter Joint

% close all;
more on;
plotResult=0;
do_setup;
tic;
cycle=1;
maxCycle=1;
% fid = fopen('betaPareto_Tensor.txt','a');
x_level=cell(length(N),1);
while(cycle<=maxCycle)

for it=8;%-6:4
Beta=10^(it);    

% for level=1:length(N)-1
%     current_n=N(level);
%     maxiter=1;
%     if(cycle>1 & current_n==N(1))
%     x0=xstar;
%     end
%     fprintf('===================== %d\n',current_n)
%     optXTM_XRF_Tensor;
%     xs{level}=xstar;
%     g_level{level}=g;
%     if(current_n~=N(end))
%         x0=downdate(xstar,1);
%     end
% end
for level=length(N):-1:1
    current_n=N(level);
    if(current_n==N(end))
        if(length(N)==1)
            maxiter=7;
        else
        maxiter=20;
        end
    else
        disp('====================== Start post-smoothing')
        maxiter=1;
        if(current_n==N(1))
            maxiter=5;
        end
    end
    if(current_n==N(1))
        plotResult=1;
    else
        plotResult=0;
    end
    fprintf('===================== %d\n',current_n)
    optXTM_XRF_Tensor;
    x_level{level}=xstar;
    save x_level x_level;
    if(current_n~=N(1))
        e=update(xstar-x0,1);
        descent1=g_level{level-1}'*e;
        descent2=g_level{level-1}'*(update(xstar)-xs{level-1});
        if(descent1<=0)
            x0=xs{level-1}+e;
        elseif(descent2<=0)
            x0=update(xstar,1);
        else
            disp('No Descent Direction');
            x0=xs{level-1};
        end
    end
end
cycle=cycle+1;
err_h=norm(xstar-W(:));
T=toc;
report_results(N);
if(Joint==1)
[f,g,f1,f2]=feval(fctn,xstar);
% fprintf(fid,'%12.4e    %12.4e     %12.4e    %12.4e      %12.4e     %f\n',Beta,f,f1,f2,err_h,T);
end

end

end
%   fclose(fid);
