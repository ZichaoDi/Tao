global dimActive
global n_delta maxiter W0 sinoS
global x_iter fiter nit

synthetic=1;
for thetaTest=1;
    tol=1e-10;
    if(thetaTest==1)
        clear xs;
        the=[20 10 20 25 ];%30 40];% 50 60];
        N=50;
        for itest=1%:length(the)
            numThetan=the(itest);
            ex_cor;
            sinoS=squeeze(-log(DisR(:,:,1)./I0'));
            b=sinoS(:);
            % null=numThetan+N^2-numThetan*(nTau+1);
            ind=find(W(:)>tol);
            ind_nempty=find(sum(L_cr(:,:,1),2)>0);
            Nzero=length(find(W(:)<=tol));
            %%============================
            Lmap=L_cr(:,:,1);
            ind_d=find(sinoS(:)<tol);
            ind_x=[];
            for i=1:length(ind_d)
                ind_x=[ind_x,find(Lmap(ind_d(i),:)>tol)];
            end
            ind_x=unique(ind_x);
            ind_nd=setdiff(1:size(Lmap,1),ind_d)';
            ind_nx=setdiff(1:N^2,ind_x)';
            l=Lmap(ind_nd,ind_nx);
            z=null(l);
            % rank0=rank(full(Lmap'*Lmap));
            % rank1=rank(full(l'*l));
            % rankCount(itest,:)=[numThetan*(nTau+1),N^2,rank0,size(l,1),size(l,2),rank1];
            %%============================
            NF = [0*N; 0*N; 0*N];
            maxiter=500;
            Lmap=sparse(squeeze(L_cr(:,:,1)));
            bounds=1;
            fctn=@(x)sfun_radon(x,sinoS,sparse(Lmap));% on attenuation coefficients miu;
            W0=W(:);
            err0=norm(W0);
            ad=linspace(-100,100,5000);
            for it_inner=1:length(ad)
                x0=0*W(:);
                % x0(ind_nx)=W(ind_nx)+ad(it_inner)*rand(size(sum(z,2)));
                x0(ind_nx)=W(ind_nx)+ad(it_inner)*sum(z,2);
                [f0,g0]=feval(fctn,x0);
                % f1=feval(fctn,x1);
                ft(it_inner,:)=[f0 norm(g0)];
                % x0(ind_nx)=l'*((l*l')\b(ind_nd));

                %% ============== optimization
                % low=zeros(N^2,1);
                % up=inf*ones(N^2,1);
                % [x,f,g,ierror] = tnbc(x0,fctn,low,up); % algo='TNbc';
                % err0=norm(x-W0,1);
                % [x_nb,f_nb,g] = tn(x0,fctn); % algo='TNbc';
                % err1=norm(x_nb-W0,1);
                % xs(:,itest)=x;
                % xs_nb(:,itest)=x_nb;
                % ft(itest,it_inner,1:5)=[f f_nb f0 err0 err1];
                %% ===============================
            end
        end

        % save('pert_initial.mat','ad','ft','xs','xs_nb');
        % save('nullInves_rank.mat','rankCount');

    elseif(thetaTest==0)
        the=[20 30 40 50];
        numThetan=10;
        for itest=1:length(the)
            N=the(itest);
            ex_cor;
            % null=numThetan+N^2-numThetan*(nTau+1);
            ind=find(W(:)>tol);
            ind_nempty=find(sum(L_cr(:,:,1),2)>0);
            Nzero=length(find(W(:)<=tol));
            opt_dr;
            err=ssim(W,align(xCOR(numThetan+1:end),W));
            tableCount(itest,:)=[null, Nzero,length(find(xCOR(numThetan+1:end)<=tol)),dimActive,err,f];
            xs{itest}=xCOR(numThetan+1:end);
            w{itest}=W;
            rankCount(itest,:)=[numThetan*(nTau+1),N^2,length(ind_nempty),length(ind),err];
            [~,stemp,v]=svd(L_cr(ind_nempty,ind,1));
            s{itest}=diag(stemp);
        end
        % save(['nullInves_ntau',num2str(numThetan),'.mat'],'tableCount','xs','w');
        save(['nullInves_ntau_rank',num2str(numThetan),'.mat'],'rankCount','s');
    end
end
% legend('d(null)','# of 0 in x^*','# of 0 in x^k','# of active set','SSIM','f')
