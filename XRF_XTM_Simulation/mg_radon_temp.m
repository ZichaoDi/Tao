global WS level
L=L_level{1};
W=W_level{1};
WS=W;
numThetan=numThetan_level(1);
nTau=nTau_level(1);
b{1}=-log(reshape(DisR',numThetan*(nTau+1),1)./I0);
x0=zeros(N(1)^2,1);
cycle=1;
maxOut=10;
xh=cell(n_level-1,1);
res=[];
x=x0;
% solver='RART';
if(strcmp(solver,'Gauss_Seidel'))
    res(1,:)=[0,norm(L'*L*x-L'*b{1}),norm(W(:)-x)];
elseif(strcmp(solver,'ART'))
    res(1,:)=[0,norm(L*x-b{1}),norm(W(:)-x)];
elseif(strcmp(solver,'RART'))
    res(1,:)=[0,norm(L*x-b{1}),norm(W(:)-x)];
end

k1=2;
k2=1e3;
while(cycle<=maxOut);
    disp('===================== Pre-smoothing')
    for level=1:n_level-1
        L=L_level{level};
        if(level>1)
            x=zeros(size(L,2),1);
        end
        if(strcmp(solver,'Gauss_Seidel'))
            xh{level}=Gauss_Seidel(L'*L,x,L'*b{level},k1);
        elseif(strcmp(solver,'ART'))
            xh{level}=solver_ART(L,x,b{level},k1);
        elseif(strcmp(solver,'RART'))
            opts=[k1,1,1e-8,0];
            xh{level}=kaczmarz(L,b{level},x,opts);
        end
        rh=b{level}-L*xh{level};
        rH=downdate_radon(rh,1,numThetan_level(level),nTau_level(level));
        b{level+1}=rH;
    end
    disp('===================== Coarsest Level Recursion')
    LH=L_level{end};
    if(strcmp(solver,'Gauss_Seidel'))
        eH=Gauss_Seidel(LH'*LH,zeros(size(LH,2),1),LH'*rH,k2);
    elseif(strcmp(solver,'ART'))
        eH=solver_ART(LH,zeros(size(LH,2),1),rH,k2);
    elseif(strcmp(solver,'RART'))
        opts=[k2,1,1e-8,0];
        eH=kaczmarz(LH,rH,zeros(size(LH,2),1),opts);
    end
    eh=eH;
    for level=n_level-1:-1:1
        L = L_level{level};
        % eh=update(eH,1);
        disp('===================== Post-smoothing')
        if(level>1)
            [eh,alpha]=lin_linear(update(eh,1),xh{level},L,b{level});
        else
            eh=update(eh,1);
            [x,alpha]=lin_linear(eh,xh{level},L,b{level});
            if(strcmp(solver,'Gauss_Seidel'))
                [x,iter]=Gauss_Seidel(L'*L,x,L'*b{level},k1);
            elseif(strcmp(solver,'ART'))
                [x,iter]=solver_ART(L,x,b{level},k1);
            elseif(strcmp(solver,'RART'))
                opts=[k1,1,1e-8,0];
                [x,iter]=kaczmarz(L,b{level},x,opts);
            end
        end
        if(level==1)
            figure(10);
            subplot(2,2,1),imagesc(reshape(eh,N(1),N(1))); colorbar;
            title(['cycle #',num2str(cycle),'updated direction']);
            subplot(2,2,2);imagesc(reshape(W(:)-xh{1},N(1),N(1))); 
            title('step to solution');colorbar;
            subplot(2,2,3);plot(res(:,1),res(:,2),'r.-');drawnow;
            e_star=W(:)-xh{1};
            subplot(2,2,4);plot(1:N(1)^2,eh,'r.-',1:N(1)^2,W(:)-xh{1},'b.-');legend('updated direction','optimal direction');drawnow;
        end
    end
    res(cycle+1,:)=[(sum(2*k1./(4.^[0:n_level-1]))+k2/24)*cycle,iter(end,2),iter(end,3)];
    cycle=cycle+1;
end
CR_mg=(res(end,2)/res(1,2))^(1/cycle);
fprintf('Convergence Rate for MG is %e \n',CR_mg);
if(strcmp(solver,'RART'))
    opts=[1e2,1,1e-8,0];
    [x,sing_rec,CR]=kaczmarz(L_level{1},b{1},x0,opts);
elseif(strcmp(solver,'ART'))
    [x,sing_rec,CR] = solver_ART(L_level{1},x0,b{1},1e4); 
elseif(strcmp(solver,'Gauss_Seidel'))
    [x,sing_rec,CR]=Gauss_Seidel(L_level{1}'*L_level{1},x0,L_level{1}'*b{1},1e2);
end
fprintf('Convergence rate for single level is %e \n',CR);
figure,plot(res(:,1),res(:,2),'r.-',sing_rec(:,1),sing_rec(:,2),'b.-');legend('mg','single-level');
