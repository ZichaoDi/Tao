%%% simple multigrid
global x0 current_n N NF
NumElement=2;
plotResult=0;
N=[17];% 9 5 3];
NF = [0*N; 0*N; 0*N];
% WH=ones(N(end),N(end),NumElement);
% for tsub=1:NumElement
% WH(:,:,tsub)=tsub*2e-1;
% end
% W=WH;

Wh=ones(N(1),N(1),NumElement);
for tsub=1:NumElement
    Wh(:,:,tsub)=tsub*2e-1;
end
W=Wh;
W_level=cell(length(N),1);
W_level{1}=W;
rng('default');
x0=1*10^(-1)*rand(size(W(:))); % Initial guess for W
xinitial=x0;
tic;
cycle=1;
maxCycle=1;
while(cycle<=maxCycle)
for level=1:length(N)-1
    current_n=N(level);
    maxiter=1;
    if(cycle>1 & current_n==N(1))
    x0=xstar;
    end
    optXTM_XRF_Tensor;
    if(current_n~=N(end))
        x0=downdate(xstar,1);
        if(cycle==1)
        W_level{level+1}=reshape(downdate(W(:),1),N(level+1),N(level+1),NumElement);
        end
        W=W_level{level+1};
    end
end


for level=length(N):-1:1
    current_n=N(level);
    if(current_n==N(end))
        if(length(N)==1)
            maxiter=1000;
        else
        maxiter=10;
        end
    else
        disp('start post-smoothing')
        maxiter=1;
    end
    if(current_n==N(1))
        plotResult=1;
    else
        plotResult=0;
    end
    optXTM_XRF_Tensor;
    if(current_n~=N(1))
        x0=update(xstar,1);
        W=W_level{level-1};%reshape(update(W(:),1),N(level-1),N(level-1),NumElement);
    end
end
cycle=cycle+1;
end
toc;
err_h=norm(xstar-W(:))/norm(xinitial-W(:))

report_results(N);