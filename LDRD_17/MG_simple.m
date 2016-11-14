%%% simple multigrid
global x0 current_n N;
global maxiter Joint

% close all;
more on;
plotResult=0;
do_setup_mg;
tic;
cycle=1;
maxCycle=2;
while(cycle<=maxCycle)

    for level=1:length(N)
        current_n=N(level);
        maxiter=1;
        fprintf('===================== %d\n',current_n)
        disp('====================== Pre-smoothing')
        if(level~=N(end))
        fctn=@(W)sfun_radon(W,xtm_level{level},I0,L_level{level});
        low=0*ones(current_n^2,1);
        up=inf*ones(current_n^2,1);
        if(current_n==N(end))
            disp('====================== Recursion')
            maxiter=20;
        else
            disp('====================== Post-smoothing')
            maxiter=1;
        end
        if(bounds)
            [x0,f,g,ierror] = tnbc (x0,fctn,low,up);
        else
            [x0,f,g,ierror] = tn (x0,fctn);
        end
        if(current_n~=N(end))
            x0=downdate(x0,1);
        end
    end
    for level=length(N):-1:1
        if(current_n~=N(1))
            e=update(xstar-x0,1);
            descent1=g_level{level}'*e;
            descent2=g_level{level}'*(update(xstar)-xs{level-1});
            if(descent1<=0)
                x0=xstar+e;
            elseif(descent2<=0)
                x0=update(xstar,1);
            else
                disp('No Descent Direction');
                x0=xstar;
            end
        end
        if(current_n==N(1))
            plotResult=1;
        else
            plotResult=0;
        end
        fprintf('===================== %d\n',current_n)
        fctn=@(W)sfun_radon(W,xtm_level{level},I0,L_level{level});
        low=0*ones(current_n^2,1);
        up=inf*ones(current_n^2,1);
        if(bounds)
            [xstar,f,g,ierror] = tnbc (x0,fctn,low,up);
        else
            [xstar,f,g,ierror] = tn (x0,fctn);
        end
        current_n=N(level-1);
    end
cycle=cycle+1;
end
err_h=norm(xstar-W(:));
T=toc;
report_results(N);
