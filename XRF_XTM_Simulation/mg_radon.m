global WS level
L=L_level{1};
W=W_level{1};
numThetan=numThetan_level(1);
nTau=nTau_level(1);
WS=W;
b=-log(reshape(DisR',numThetan*(nTau+1),1)./I0);
x0=zeros(N(1)^2,1);
cycle=1;
maxOut=10;

res=[];
x=x0;
% solver='RART';
if(strcmp(solver,'Gauss_Seidel'))
    res(1,:)=[0,norm(L'*L*x-L'*b),norm(W(:)-x)];
elseif(strcmp(solver,'ART'))
    res(1,:)=[0,norm(L*x-b),norm(W(:)-x)];
elseif(strcmp(solver,'RART'))
    res(1,:)=[0,norm(L*x-b),norm(W(:)-x)];
end

k1=2;
k2=1e3;
while(cycle<=maxOut);
    disp('===================== Pre-smoothing')
    if(strcmp(solver,'Gauss_Seidel'))
        xh=Gauss_Seidel(L'*L,x,L'*b,k1);
    elseif(strcmp(solver,'ART'))
        xh=solver_ART(L,x,b,k1);
    elseif(strcmp(solver,'RART'))
        opts=[k1,1,1e-8,0];
        xh=kaczmarz(L,b,x,opts);
    end
    rh=b-L*xh;
    rH=downdate_radon(rh,1);
    disp('===================== Recursion')
    if(strcmp(solver,'Gauss_Seidel'))
        eH=Gauss_Seidel(LH'*LH,zeros(size(LH,2),1),LH'*rH,k2);
    elseif(strcmp(solver,'ART'))
        eH=solver_ART(LH,zeros(size(LH,2),1),rH,k2);
    elseif(strcmp(solver,'RART'))
        opts=[k2,1,1e-8,0];
        eH=kaczmarz(LH,rH,zeros(size(LH,2),1),opts);
    end
    eh=update(eH,1);
    [x,alpha]=lin_linear(eh,xh,L,b);
    disp('===================== Post-smoothing')
    if(strcmp(solver,'Gauss_Seidel'))
        [x,iter]=Gauss_Seidel(L'*L,x,L'*b,k1);
    elseif(strcmp(solver,'ART'))
        [x,iter]=solver_ART(L,x,b,k1);
    elseif(strcmp(solver,'RART'))
        opts=[k1,1,1e-8,0];
        [x,iter]=kaczmarz(L,b,x,opts);
    end
    res(cycle+1,:)=[(2*k1+k2/24)*cycle,iter(end,2),iter(end,3)];
    figure(10);
    subplot(2,2,1),imagesc(reshape(eh,N(1),N(1))); colorbar;
    title(['cycle #',num2str(cycle),'updated direction']);
    subplot(2,2,2);imagesc(reshape(W(:)-xh,N(1),N(1))); 
    title('step to solution');colorbar;
    subplot(2,2,3);plot(res(:,1),res(:,2),'r.-');drawnow;
    e_star=W(:)-xh;
    subplot(2,2,4);plot(1:N(1)^2,map1D(alpha*eh,e_star),'r.-',1:N(1)^2,W(:)-xh,'b.-');drawnow;
    cycle=cycle+1;
end
CR_mg=(res(end,2)/res(1,2))^(1/cycle);
fprintf('Convergence Rate for MG is %e \n',CR_mg);
if(strcmp(solver,'RART'))
    opts=[1e2,1,1e-8,0];
    [x,sing_rec,CR]=kaczmarz(L,b,x0,opts);
elseif(strcmp(solver,'ART'))
    [x,sing_rec,CR] = solver_ART(L,x0,b,1e4); 
elseif(strcmp(solver,'Gauss_Seidel'))
    [x,sing_rec,CR]=Gauss_Seidel(L'*L,x0,L'*b,1e2);
end
fprintf('Convergence rate for single level is %e \n',CR);
figure,plot(res(:,1),res(:,2),'r.-',sing_rec(:,1),sing_rec(:,2),'b.-');legend('mg','single-level');
