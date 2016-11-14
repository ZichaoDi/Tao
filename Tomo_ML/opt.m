fctn=@(MU)sfun_radon(MU,Mt,I0,L);% on attenuation coefficients miu;

rng('default');
x0 = 0*10^(0)*rand(m(1),m(2),NumElement);
x0=x0(:);
err0=norm(W0-x0);
NF = [0*N; 0*N; 0*N];
maxOut=100;

