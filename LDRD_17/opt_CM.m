function [x] = opt_CM(DisR_Simulated,I0)
global maxiter NF N 
Mt=-log(DisR_Simulated/I0);
fctn=@(x)CenterMass(x,Mt);
NF = [0*N; 0*N; 0*N];
rng('default');
x0 = 0*rand(3,1);
x=lsqnonlin(fctn,x0);
maxiter=100;
[x_tn,f,g,ierror] = tn (x0,fctn);


function [f,g]=CenterMass(c,DisR);
global thetan numThetan nTau
m0=sum(DisR(:))/(numThetan*(nTau+1));
for n=1:numThetan
    tc(n)=sum(DisR(:,n).*[1:nTau+1]')/m0;
end
theta=thetan/180*pi;
tcp=c(1)*cos(theta)+c(2)*sin(theta)+c(3);
f=sum((tcp-tc).^2)/2;
% figure, plot(1:numThetan,tcp,'r.-',1:numThetan,tc,'b.-')
g=[sum((tcp-tc).*cos(theta));sum((tcp-tc).*sin(theta));sum(tcp-tc)];
