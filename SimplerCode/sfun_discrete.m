function [f,g]=sfun_discrete(I,M,I0,Ltol,thetan,omega,m,dTau)
% load DisCornerSphere20theta90;
%%===== Reconstruction discrete objective 
%%===== I: current intensity approximation
%%===== L: intersection length matrix
%%===== R: Radon transform with t beam lines and theta angles
%%===== f: sum_i ||e^T(L_i.*I)e-Ri||^2, i=1..theta
I=reshape(I,size(I0));
e=ones(size(I,1),1);
nTau=ceil(sqrt(omega(2)^2+omega(4)^2)/dTau);
f=0;
g=zeros(m(1),m(2));
J=g;
for n=1:length(thetan)
    theta=thetan(n)/180*pi;
    if(theta<pi/2)
Mt=flipud(M(:,n));
    else
        Mt=M(:,n);
    end
    if(theta==0 | theta==pi)
        diagL1=ceil(omega(1));
        diagL2=ceil(omega(2));
    elseif(theta==pi/2)
        diagL1=ceil(omega(3));
        diagL2=ceil(omega(4));
    elseif(theta<pi/2)
        diagL1=ceil(omega(2)*tan(theta)/dTau);
        diagL2=ceil(omega(4)/dTau);
    else
        diagL1=0;%ceil(omega(2)*tan(theta)/dTau);
        diagL2=ceil((omega(4)-omega(2)*tan(theta))/dTau);
    end
    Ttol=linspace(-diagL1+0,diagL2-0,nTau);
    sum_Tau=0;
    for i=1:length(Ttol)
        
        L=Ltol{n,i};
        if(isempty(L))
            L=zeros(m(1),m(2));
        end
        Rdis=e'*(I.*L)*e; %% Discrete case
        sum_Tau=sum_Tau+(Rdis-Mt(i))^2;
        g=g+2*(Rdis-Mt(i)).*L;
        J=J+L;
    end
    f=f+sum_Tau;
end
condN=cond(J);
[U,S,V]=svd(J);
minS=min(diag(S));
maxS=max(diag(S));
fprintf('Condition %d, maximumS %d, minimumS %d \n',condN,maxS, minS);
g=g(:);