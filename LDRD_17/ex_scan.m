global initialize
global NF W0 err0 maxiter
maxiter=50;

Ntest=2;
DisR=[];
L_cr=[];
nslice=150;
initialize=1;
ind_cr=1;
n_delta=0;
for ind_scan=1:Ntest;
        do_setup; 
        L_cr(:,:,ind_scan)=L;
        DisR(:,:,ind_scan)=DisR_Simulated;
end
NF = [0*N; 0*N; 0*N];
Mt=-log(DisR./I0);
%% make sure sample area is positive
touchedInd1=zeros(numThetan,nTau+1,N^2);
touchedInd2=zeros(numThetan,nTau+1,N^2);
for n=1:numThetan
    for i=1:nTau+1
        % MU_pert=map1D(w_drift(:,:,i),[0,1]);
        touchedInd1(n,i,find(W(:)>0&L_cr((nTau+1)*(n-1)+i,:,1)'>0))=1;
        touchedInd2(n,i,find(W(:)>0&L_cr((nTau+1)*(n-1)+i,:,2)'>0))=1;
    end
end
touched1=reshape(sum(sum(touchedInd1,1),2),N,N);
touched2=reshape(sum(sum(touchedInd2,1),2),N,N);
% x1=L_cr(:,:,1)\reshape(Mt(:,:,1)',numThetan*(nTau+1),1);
% x2=L_cr(:,:,2)\reshape(Mt(:,:,2)',numThetan*(nTau+1),1);
L_cr=L_cr./scale;
W0=W(:);
err0=norm(W0);
fctn1=@(x)sfun_radon(x,Mt(:,:,1)',L_cr(:,:,1));% on attenuation coefficients miu;
[x1,f,g,ierror] = tnbc (ones(N^2,1),fctn1,zeros(N^2,1),inf*ones(N^2,1)); % algo='TNbc';
fctn2=@(x)sfun_radon(x,Mt(:,:,2)',L_cr(:,:,2));% on attenuation coefficients miu;
[x2,f,g,ierror] = tnbc (zeros(N^2,1),fctn2,zeros(N^2,1),inf*ones(N^2,1)); % algo='TNbc';
nrow=3;
clims=[0,1];
figure,subplot(nrow,2,1), imagesc(touched1);colorbar; 
subplot(nrow,2,2),imagesc(touched2);colorbar; 
subplot(nrow,2,[3 4]),plot(1:N,touched1(:,round(N/2)),'r.-',1:N,touched2(:,round(N/2)),'b.-'); legend('error-free','perturbed');
subplot(nrow,2,5);imagesc(reshape(x1,N,N));colorbar; 
subplot(nrow,2,6);imagesc(reshape(x2,N,N));colorbar;


