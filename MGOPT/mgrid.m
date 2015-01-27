function [v, varargout] = mgrid(v0,fnl,res_prob,step_bnd)
%--------------------------------------------------
% Multigrid algorithm (McCormick, pp. 66-67)
%--------------------------------------------------
% Usage: v = mgrid(v0,fnl,res_prob)
%--------------------------------------------------
global current_fnl
global N current_n
global bounds v_low v_up
global GRAPH_N_OLD GRAPH_INDEX_OLD it L lambdatn lambdaind comp shiftc
global Hess_norm_H % norm of coarse-grid Hessian 
global test_v_in_fmincon v_after_r row bind spec_bd IHh IhH vvl coalow vvu bindH
row=1;
%----------------------------------------------------------------------
% DECIDE WHETHER TO PRINT ASSESSMENT TEST RESULTS
% AND WHETHER TO PAUSE AFTER PRINTING TESTS
%----------------------------------------------------------------------
assess_print = 0; % 0=false, 1=true
assess_pause = 0; % 0=false, 1=true
%----------------------------------------------------------------------
% UPDATE MULTIGRID GRAPH
%----------------------------------------------------------------------
h_MG_graph = findobj('Tag', 'multigrid_graph');
figure(h_MG_graph);

GRAPH_N_NEW = current_n;
IND_OLD = find(N==GRAPH_N_OLD);
IND_NEW = find(N==GRAPH_N_NEW);
plot([GRAPH_INDEX_OLD GRAPH_INDEX_OLD+1],[IND_OLD IND_NEW]);
plot([GRAPH_INDEX_OLD GRAPH_INDEX_OLD+1],[IND_OLD IND_NEW],'x');
GRAPH_N_OLD = current_n;
GRAPH_INDEX_OLD = GRAPH_INDEX_OLD+1;
%%%%%%%%%%############################ Test the variable after recursion
test_v_in_fmincon=0;
%%%%%%%%%%%%######################################################
%--------------------------------------------------
n = current_n;
global_setup(n);
alpha = 0; % added to try to deal with error message
%--------------------------------------------------
% Determine update/downdate constant
if (n>min(N));
    v_test_h = 0*v0;
    v_test_h(1) = 1;
    v_D1     = downdate(v_test_h,1); % downdated test vector (first column of D)
    v_test_H = 0*v_D1;
    v_test_H(1) = 1;
    v_U1     = update(v_test_H,1);   % updated test vector (first column of U)
    IhH_con = v_U1(1)/v_D1(1);
    if (assess_print)
        fprintf('--------------------------------------------------------------\n');
        fprintf('Update/Downdate Constant: %d\n', IhH_con);
        fprintf('--------------------------------------------------------------\n');
        if (assess_pause); disp('Hit any key to continue'); pause; end;
    end;
end;
%--------------------------------------------------
current_fnl = fnl;
j           = find(N==current_n);
nmin        = N(end);
if (bounds);
    if(res_prob & spec_bd)
        [v_low,v_up] = set_bounds2(j,bindH);
        vvl{j}=v_low;
        vvu{j}=v_up;
        
    else
        [v_low,v_up] = set_bounds(j,res_prob);
        vvl{j}= v_low;
        vvu{j}=v_up;
    end
    %[v_low,v_up] = set_bounds(j);
end;
%--------------------------------------------------
if (res_prob);
    myfun = 'sfun_mg';
else
    myfun = 'sfun';
end;
%--------------------------------------------------
T  = ['In  mgrid: n = ' num2str(current_n)];
disp(T)
%--------------------------------------------------
if (res_prob);
    [F_on_entry, G_on_entry] = sfun_mg(v0);
else
    [F_on_entry, G_on_entry] = sfun(v0);
end;
npts = numel(v0) - 1;
%--------------------------------------------------
if (current_n <= nmin);
    %--------------------------------------------------
    % solve (exactly) problem on coarsest grid
    %--------------------------------------------------
    nit_solve = 250;
    if (step_bnd == 0);
        if (bounds);
            
            nit = nit_solve; [v,F,G,ierror]  = tnbcm(v0,myfun,v_low,v_up,nit);

            coalow=v_low;
            
        else
            nit = nit_solve; [v,F,G,ierror]  = tnm  (v0,myfun,nit);
        end;
    else
        [low1,up1] = get_bnd (v0,step_bnd);
        nit = nit_solve; [v,F,G,ierror]  = tnbcm(v0,myfun,low1,up1,nit);
        end
%                               load vs17l01; vc=vs17l01;
%                             figure(133);
%                             plot(1:current_n,v,'r.-',1:current_n,v0,'b.-',1:current_n,vc,'k.-'); title('red:vs; blue: v0; black: vc;')
            %              lamH=zeros(size(v0));
            % %             %%%=========================================
            %                                     indl=find(v<=v_low);
            %                                     lamH(indl)=G(indl);
            %                                     load lam1l257;
            %                                     lams=lam1l257;
            % %             %%%=========================================
            % %             indu=find(v>=v_up);
            % %             lamH(indu)=-G(indu);
            % %             load lam129us; lams=lam129us;
            % %             load lamH65;
            % %             %%%=========================================
            %              lamh = update (lamH, 1);
            %             figure(23);
            % %             % subplot(2,1,1);
            %            plot(1:N(1),lams,'ro-',1:N(1),lamh,'b*-'); title('red: fine lambda; blue: updated coarse lam')
            %             %   subplot(2,1,2);plot(1:N(1),lamh,'b*-');
            %             %plot(1:N(2),lamH,'b.-');%,1:N(2),lamH65,'ro-');%             title('red: original coarse lam; blue: shifted coarse lam')
            %%%%%%%%%%%%%%%%###########################################################
            %          figure(4444);
            %             subplot(2,1,1)
            % %             load vs1s7
            % %             load vs1s
            %             plot(1:N(j),v_low,'b.-');title(' pink: v-current_v')
            %             subplot(2,1,2)
            %             plot(1:N(j-1)^2,vs1s7-bind{j},'g.-');title(' green: vsh-vh')
            %             figure(12);
            %             subplot(2,1,1);
            %             plot(1:N(j)^2,bindH{j},'g.-');title('blue:v; green: current_v; pink: vsH; red: v_low;')
            %             subplot(2,1,2);
            %             plot(1:N(j-1)^2,bind{j},'r.-'); title('red: vh; blue: vsh; green: v_low_h')
            %             pause
            %%%%%%%%%%%%%%%%###############################################
    
else
    %--------------------------------------------------
    % "relax" on current grid
    %--------------------------------------------------
    if (step_bnd == 0);
        if (bounds);
            nit = 1; [v,F,G,ierror,eig_val]  = tnbcm(v0,myfun,v_low,v_up,nit);
        else
            nit = 1; [v,F,G,ierror,eig_val]  = tnm  (v0,myfun,nit);
        end;
    else
        [low1,up1] = get_bnd (v0,step_bnd);
        nit = 1; [v,F,G,ierror,eig_val]  = tnbcm(v0,myfun,low1,up1,nit);
    end;
    if (res_prob);
        [F_before_recursion, G_before_recursion] = sfun_mg(v);
        [F_foo, G_before_recursion_reg] = sfun(v);
    else
        [F_before_recursion, G_before_recursion] = sfun(v);
        G_before_recursion_reg = G_before_recursion;
    end;
    %--------------------------------------------------
    % Form new "rhs" vector
    %--------------------------------------------------
    current_v   = downdate(v,3);

    if (bounds);
        if(res_prob & spec_bd)
            [v_low_1,v_up_1] = set_bounds2(j+1,bindH);
        else
            [v_low_1,v_up_1] = set_bounds(j+1,res_prob);
        end
        current_v = bound_project(current_v,v_low_1,v_up_1);
    end;
    if(bounds)
        bind{j+1}=v;
        bindH{j+1}=current_v;
        
    end;
%    figure(112);plot(2:2:N(1),current_v,'r.-',1:N(1),v,'b.-')
    [Fv  ,Gv  ] = sfun   (v);
    [Fvmg,Gvmg] = sfun_mg(v);
    %%%%#####################################
    if(bounds)
        lamG=Gv;
        lambda_low=zeros(length(v),1);
        lambda_low(lambdaind.low)=lamG(lambdaind.low);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lambda_upper=zeros(length(v),1);
        lambda_upper(lambdaind.up)=-lamG(lambdaind.up);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lambdatn=struct('lower',lambda_low,'upper',lambda_upper);
    end;
    
    %%%%#######################################
    %figure(121);plot(lambdatn.lower)
    ddGv= downdate(Gv,1);
    if(bounds==1)
        Gv=Gv-lambdatn.lower+lambdatn.upper;
    end
    dGv         = downdate(Gv,1);
    if(bounds)
        current_lambdatn_lower=downdate(lambdatn.lower,1);
        current_lambdatn_upper=downdate(lambdatn.upper,1);
        cc= current_lambdatn_upper;
        current_lambdatn=struct('lower',current_lambdatn_lower,'upper',current_lambdatn_upper);
    end
    fnl2        = downdate(fnl,1);
    %%----------------------------------------------------
    test_sep = true;
    j = find(N==current_n);
    if (test_sep)
        %-------------------------------------------------
        % Generate "probe" vector p_h
        %-------------------------------------------------
        % rml randn('state',0);
        % rml p_0h = randn(size(v));
        
        p_0h = 2*(rand(size(v)) - 0.5);
        p_iter_max = 1; % I used "5" instead of "1" here
        p_1h = p_0h;
        for jj = 1:p_iter_max;   % Iterate to accentuate the low "frequencies"
            p_1h = update(downdate(p_1h,0),0);
        end;
        p_h  = p_0h - p_1h;    % This should leave mostly high-frequencies in p_h.
        %------------------------------------------------
        % Downdate the probe vector and compute
        % the fine-grid matrix-vector product Ap_h0.
        %------------------------------------------------
        h_diff = sqrt(eps)*norm(p_h);
        [Fvh, Gvh]  = sfun(v + h_diff*p_h);
        Ap_h0 = (Gvh - Gv) / h_diff;
        %------------------------------------------------
        Ap_H = downdate(Ap_h0,0);
        Ap_h = Ap_h0; % RML hack.
        n_Ap_h = norm(Ap_h);
        %-----------------------------------------------
        % Iterate to emphasize the low frequencies.
        % If the Hessian is not separable, then the
        % low frequencies will be large; otherwise
        % they will be small.
        %-----------------------------------------------
        A_iter_max = 1; % RML used "5" here
        for jj = 1:A_iter_max
            Ap_h = update(downdate(Ap_h,0),0);
            Ap_H = downdate(update(Ap_H,0),0);
        end
        %------------------------------------------------
        % Plot the various results
        % In the general case the FFT may not be an
        % appropriate tool to analyze the results.
        %------------------------------------------------
        kkk = numel(eig_val);
        while (kkk > 0 && isempty(eig_val{kkk}))
            kkk = kkk-1;
        end
        if (kkk > 0)
            Hessian_norm = max(abs(eig_val{kkk}));
            sep_denom = norm(p_h) * Hessian_norm;
        else
            sep_denom = n_Ap_h;
        end;
    end;
    %%----------------------------------------------------
    
    j           = j+1;
    current_n   = N(j);
    %--------------------------------------------------
    n = current_n;
    global_setup(n);
    %--------------------------------------------------
    
    
    [Fv2 ,Gv2 ] = sfun   (current_v);
    %figure(25);plot(1:N(2),ddGv,'r.-',1:N(2),Gv2,'b.-')
    %     HL=find(current_v>=v_low_1);
    %     HU=find(current_v<=v_up_1);
    %       load lamH65;
    %     current_lambdatn.lower=lamH65;%0.*current_v;
    %     current_lambdatn.upper=0.*current_v;
    %     current_lambdatn.lower(HL)=Gv2(HL);
    %     current_lambdatn.upper(HU)=-Gv2(HU);
    
    % figure(24);plot(1:N(2),cc,'r.-',1:N(2),lamH65,'b.-'); title('red: coarsed lam; blue: exact coarse lam')
    %%----------------------------------------------------
    %%%#########################shifted with Lagrangian multiplier
    if(bounds)
        Gv2=Gv2-current_lambdatn.lower+current_lambdatn.upper;
    end
    tau         = Gv2 - dGv;
    shiftc1=norm(Gv2);
    shiftc2=norm(dGv);
    shiftc3=norm(tau);
    shiftc=struct('a',shiftc1,'b',shiftc2,'c',shiftc3);
    fnl2        = fnl2 + tau;
    %shiftnorm=norm(fnl2)
    %   figure(121);plot(1:N(2),Gv2,'r.-',1:N(2),dGv,'b.-'); title('shape of the shifted part')
    % figure(11);plot(fnl2)
    %--------------------------------------------------
    % Bound step on next-coarser grid
    %--------------------------------------------------
    bnd_1     = norm(fnl2,'inf');
    bnd_2     = norm( Gv2,'inf');
    bnd_3     = norm( dGv,'inf');
    step_bnd1 = 10 * max([bnd_1 bnd_2 bnd_3]);
    %########################################
    % step_bnd1;
    %     step_bnd1 = input('Enter value for step_bnd1:  ')
    no_bounds = 1;  % 1: no extra bound constraint; 0: extra bound constraint
    if (no_bounds);
        step_bnd1 = 0;
        disp(' ')
        disp('#################################')
        disp('### NO EXTRA BOUND CONSTRAINT ###')
        disp('#################################')
        disp(' ')
    end;
    %########################################
    %--------------------------------------------------
    % update solution on current grid (V cycle)
    %--------------------------------------------------
    [e2, pred]  = mgrid(current_v,fnl2,1,step_bnd1);
    
    j           = j-1;
    current_n   = N(j);
    e           = update(e2-current_v,1);
    
    if (Hess_norm_H > 0)
        sep_denom = norm(p_h) * Hess_norm_H;
        sep_quot = norm(Ap_H) / sep_denom;
        if (assess_print);
            fprintf('--------------------------------------------------------------\n');
            fprintf('SEPARABILITY TEST (nh = %d, using |G_H|) = %f\n', N(j),sep_quot);
        end;
    end
    %--------------------------------------------------
    % Compute approximate Rayleigh quotient.
    %--------------------------------------------------
    kk = numel(eig_val);
    while (kk > 0 && isempty(eig_val{kk}))
        kk = kk-1;
    end
    if (kk > 0)
        Hessian_norm = max(abs(eig_val{kk}));
        Hess_norm_H = Hessian_norm;
        accrcy = 128*eps;
        vnorm  = norm(v, 'inf');
        He = gtims (e, v, G_before_recursion_reg, accrcy, vnorm, @sfun);
        approx_Rayleigh_quot = (dot(He,He)/(Hessian_norm * dot(e,He)));
        if (assess_print);
            fprintf('MESH COMPLEMENTARITY TEST (nh = %d): %f\n',N(j),approx_Rayleigh_quot);
        end;
    else
        if (assess_print);
            fprintf('\n');
            fprintf('MESH COMPLEMENTARITY TEST\n');
            fprintf('  Comparing meshes nh = %d and nH = %d.\n', N(j), N(j+1));
            fprintf('  The test could not be performed because no estimate of the norm\n');
            fprintf('  of the Hessian was available!  This is because modlnp returned\n');
            fprintf('  before the requisite Lanczos information was collected!\n');
        end;
        approx_Rayleigh_quot = -1000; % This is a dummy value to indicate test failure
        Hessian_norm = 0; % This is a dummy value to indicate test failure
        Hess_norm_H = Hessian_norm;
    end
    %----------------------------------------------------------------------
    % UPDATE MULTIGRID GRAPH
    %----------------------------------------------------------------------
    figure(findobj('Tag', 'multigrid_graph'));
    GRAPH_N_NEW = current_n;
    IND_OLD = find(N==GRAPH_N_OLD);
    IND_NEW = find(N==GRAPH_N_NEW);
    plot([GRAPH_INDEX_OLD GRAPH_INDEX_OLD+1],[IND_OLD IND_NEW]);
    plot([GRAPH_INDEX_OLD GRAPH_INDEX_OLD+1],[IND_OLD IND_NEW],'x');
    GRAPH_N_OLD = current_n;
    GRAPH_INDEX_OLD = GRAPH_INDEX_OLD+1;
    %--------------------------------------------------
    n = current_n;
    global_setup(n);
    if (bounds);
        if(res_prob & spec_bd)
            [v_low,v_up] = set_bounds2(j,bindH);
        else
            [v_low,v_up] = set_bounds(j,res_prob);
        end
    end;
    %--------------------------------------------------
    current_fnl = fnl;
    if(bounds)
        [Fvmg,Gvmg]=merit_tn(v,Fvmg,Gvmg);
        %%%===========================================================
        ee=e;
        indu=find(v>=v_up);
        for t=1:length(indu)
            if (e(indu(t))>0) e(indu(t))=0; end
        end
                indl=find(v<=v_low);
        for t=1:length(indl)
            if (e(indl(t))<0) e(indl(t))=0; end
        end
        %%%%%%========================================================
    end
    
    testd       = e'*Gvmg;
    
        
    %======================================================================
    % plot search direction, step to solution
    if (current_n==N(1) & comp)
        load vs575u15;             vstar=vs575u15;
        %         load vc
        %             load vs1s31;
        %             vsc=vs1s31;
        %                     load vs_129_u; vstar = vs_129_u;
        %         load vs_129_b; vstar = vs_129_b;
        %  aprerr=downdate(vstar-v,1);
        %             figure(13);
        %             plot(1:31^2,e2-current_v,'r.-');
        %         title('Search direction_blue, coarse step to solution_red (coarse)')
        %         subplot(5,2,it);
        %             %             plot(1:N(2),current_v+inv(IhH*IHh)*aprerr,'b.-',1:N(2),e2,'r.-',1:N(2),vsc,'g.-'); %title ('Blue: Approximation of shifted vstar by fine search direction; Red: real shifted vstar (coarse); Green: solution of original coarse problem')
        %           figure(19);
        % %             %             subplot(5,2,it);
        %  vstarH=downdate(vstar,0);
        %  subplot(2,1,1);
        %   plot(1:N(2),vstarH,'b.-');title('downdate of vstar_blue v.s. e2_red')
        %   subplot(2,1,2);plot(1:N(2),e2,'r.-');
        % %  figure(26);
        %  vstarh=update(e2,0);
        % subplot(2,1,1);
        %  plot(1:N(1)^2,vstar,'b.-');title('vstar v.s. update of e2')
        %  subplot(2,1,2);plot(1:N(1)^2,vstarh,'r.-');
        %title('Search direction: Approximation of shifted vstar-blue, real shifted vstar (coarse)-red; original coarse solution-green')
        %          figure(18);plot(e,'r.-'); title('Search direction (fine)')
        %         tt = linspace(0,1,n+2)'; tt = tt(2:end-1);
        
        %             test=(vstar-v)'*Gvmg
        figure(16);
        % subplot(5,1,it);
%         e(1:27)=ee(1:27);
%         e(558:end)=ee(558:end);
        plot(1:current_n,vstar-v,'b.-',1:current_n,e,'r.-',1:current_n,ee,'k.-'); %for i=1:254 text(tt(i),e(i),num2str(i));end ; for i=1:254 text(tt(i),vvt(i),num2str(i));end ;
        %text(58,e(58),num2str(58)); text(62,e(62),num2str(62));
        title('Step to solution (blue), Search direction (red)')
        figure(17);
        subplot(1,2,1); plot(vstar-v,'b.-'); title('step to solution')
        %plot(1:N(1),v,'b.-',1:N(1),vstar,'r.-');title('blue: current variable; red: solution');
        subplot(1,2,2);plot(e,'r.-');title('search direction')
        %         figure(20);
        %         subplot(5,1,it)
        %         plot(1:N(2),e2,'g.-',1:N(2),coalow,'r.-',1:N(2),vc{it},'b.-'); title('Current estimate of solution')
        % vc{it}=e2;
        % save vc vc;
        %         subplot(1,2,2)
        %         plot(vstar,'k'); title('Solution')
        
        %         disp('Hit any key to continue')
        %         pause(1)
    end
    %%=====================================================================
    
    if (testd<0)
        if (bounds)
            
            alpha0 = stpmax1 (v, e, v_low, v_up);
        else
            alpha0 = 1;
        end;
        ind_test_dir = 0; % indicator: if =1 then plot v and v+e
        if (ind_test_dir)
            nn = length(v);
            tt = linspace(0,1,nn);
            figure(701)
            plot(tt,v), hold on
            plot(tt,v+e,'r'), hold off
            disp('Hit any key to continue')
            pause
        end
        if(bounds==0)
            [v, Fv, Gv, nf1, alpha, ierror, dfdp] = ...
                lin2 (e, v, Fvmg, alpha0, Gvmg, myfun);
        elseif(bounds==1)   %%%%%%%%%%%%%This line search is based on merit function
            [v, Fv, Gv, nf1, alpha, ierror, dfdp] = ...
                lin3 (e, v, Fvmg, alpha0, Gvmg, myfun);
        end
        if (res_prob);
            [F_after_linesearch, G_after_linesearch] = sfun_mg(v);
        else
            [F_after_linesearch, G_after_linesearch] = sfun(v);
        end;
        dfdp = dfdp/npts;
        if (assess_print);
            fprintf('NONLINEARITY TEST (nh = %d, alpha = %8.2e): %8.2e\n',N(j),alpha,dfdp);
        end;
    else
        disp('No descent direction: no step taken')
        v  = v + 0*e;
    end;
    
    v_after_r=v;
    
    if (res_prob);
        [F_after_recursion, G_after_recursion] = sfun_mg(v);
    else
        [F_after_recursion, G_after_recursion] = sfun(v);
    end;
    
    aared = F_before_recursion - F_after_recursion;
    mct_pred   = IhH_con*pred;   % predicted reduction ("raw")
    mct_pred_a = alpha*mct_pred; % predicted reduction ("pro-rated" by alpha)
    mct_act    = aared;          % actual reduction
    mct_err    = abs(mct_act-mct_pred)  /abs(mct_act); % relative error ("raw")
    mct_err_a  = abs(mct_act-mct_pred_a)/abs(mct_act); % relative error ("pro-rated")
    if (assess_print);
        fprintf('MODEL CONSISTENCY TEST: nh = %d\n', N(j));
        fprintf('  Reduction(predicted) (raw/pro-rated): %e / %e\n', mct_pred,mct_pred_a);
        fprintf('  Reduction(actual)                   : %e\n', mct_act);
        fprintf('  Relative error (raw/pro-rated)      : %f / %f\n', mct_err,mct_err_a);
        fprintf('--------------------------------------------------------------\n');
        if (assess_pause); disp('Hit any key to continue'); pause; end;
    end;
    %--------------------------------------------------
    % "relax" on current grid
    %--------------------------------------------------
    disp(T)
    test_v_in_fmincon=1;
    if (step_bnd == 0);
        if (bounds);
            
            nit = 1; [v,F,G,ierror]  = tnbcm(v,myfun,v_low,v_up,nit);
            
        else
            nit = 1; [v,F,G,ierror]  = tnm  (v,myfun,nit);
        end;
    else
        [low1,up1] = get_bnd (v,step_bnd);
        nit = 1; [v,F,G,ierror]  = tnbcm(v,myfun,low1,up1,nit);
    end;
end

%%%%%%%%%%%%################test Lagrangian multiplier after post-smoothing
% if(bounds & current_n==N(1))
%     load lambda63
%     lambdastar=lambda63;
%     [f,g]=sfun(v);
%     lamup=zeros(N(1)^2,1);
%     lamup(lambdaind.up)=-g(lambdaind.up);
%     ttt=lamup(lambdaind.up);
%     lamlow=zeros(N(1)^2,1);
%     lamlow(lambdaind.low)=g(lambdaind.low);
% %     slack_up=lamup'*(v-v_up)
% %     slack_low=lamlow'*(v_low-v)
% %     indl=lambdaind.low
% %     indu=lambdaind.up
%     figure(123); %title('error of current lambda with exact lambda')
% %    subplot(5,1,it);
%     plot(1:N(1)^2,lamup-lambdastar,'r.');%,1:N(1)^2,lambda63,'bo');%,1:N(1),g,'m*')
% end


%%%%%%%%%%%%%%%##############################################

% Compute the actual reduction seen on the current mesh level
% to serve as the predicted reduction on the next finer mesh.
% RML: The scaling by current_n needs to be corrected.
if (res_prob);
    [F_on_exit, G_on_exit] = sfun_mg(v);
else
    [F_on_exit, G_on_exit] = sfun(v);
end;

if (current_n > nmin)
    ared = (F_on_entry - F_on_exit);
else
    ared = (F_on_entry - F_on_exit);
end

if (nargout == 2)
    varargout{1} = ared;
end