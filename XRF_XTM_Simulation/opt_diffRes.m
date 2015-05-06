global low up penalty
global W0 current_n WS
global SigMa_XTM SigMa_XRF Joint
global err0 fiter nit maxiter xinitial
global Beta TempBeta

close all;
more on;
plotResult=0;
do_setup;
cycle=1;
XRF=xrf_level{1};
DisR=xtm_level{1};
nTau=nTau_level(1);
x_level_J=cell(length(N),1);
x_level_R=cell(length(N),1);
x_level_T=cell(length(N),1);
fid = fopen('Convergence_Factor_Res.txt','a');
maxiter=300;
for level=1:length(N)
    current_n=N(level);
    W0= W_level{level};
    W0=W0(:);
    L=L_level{level};
    GlobalInd=GI_level{level};
    SelfInd=SI_level{level};
    m=m_level(level,:);
    SigMa_XTM=SigmaT{level};
    SigMa_XRF=SigmaR{level};
    penalty=0;
    fctn=@(W)sfun_Tensor_Joint(W,XRF,DisR,MU_e,M,NumElement,L,GlobalInd,SelfInd,thetan,m,nTau,I0);
        xinitial=x0;
        low=0*ones(size(x0));
        up=1e6*ones(size(x0));
    for J_i=-1:1
        fprintf('===================== %d, %d\n',current_n,J_i)

        
        if(J_i==1)
            Beta=1e8; TempBeta=1;
            t0=cputime;
            [xs,f,g,ierror] = tnbc (x0,fctn,low,up);
            xstar=xs;
            t_J=cputime-t0;
            if(level>1)
                for j_sub=level:-1:2
                    xs=update(xs,1);
                end
            else
                xs=xstar;
            end
            x_level_J{level}=xstar;
            convFac_J=(fiter(end)/fiter(1))^(1/(nit+1));
            
            errTol_J=norm(xs-WS(:))/norm(err0);
        elseif(J_i==0)
            Beta=0; TempBeta=1;
            t0=cputime;
            [xs,f,g,ierror] = tnbc (x0,fctn,low,up);
            xstar=xs;
            t_R=cputime-t0;
            if(level>1)
                for j_sub=level:-1:2
                    xs=update(xs,1);
                end
            else
                xs=xstar;
            end
            x_level_R{level}=xstar;
            convFac_R=(fiter(end)/fiter(1))^(1/(nit+1));
            errTol_R=norm(xs-WS(:))/norm(err0);
        elseif(J_i==-1)
            Beta=1; TempBeta=0;
            t0=cputime;
            [xs,f,g,ierror] = tnbc (x0,fctn,low,up);
            xstar=xs;
            t_T=cputime-t0;
            if(level>1)
                for j_sub=level:-1:2
                    xs=update(xs,1);
                end
            else
                xs=xstar;
            end
            x_level_T{level}=xstar;
            convFac_T=(fiter(end)/fiter(1))^(1/(nit+1));
            errTol_T=norm(xs-WS(:))/norm(err0);
        end
    end
    if(level~=length(N))
    x0=downdate(x0,1);
    end
    fprintf(fid,'%d    %12.4e     %12.4e    %12.4e     %12.4e    %12.4e     %12.4e     %12.4e    %12.4e     %12.4e\n', current_n, convFac_J, t_J, errTol_J, convFac_R, t_R, errTol_R, convFac_T, t_T, errTol_T);
% save('xstar_res.mat','x_level_J','x_level_R','x_level_T')
end

fclose(fid);

