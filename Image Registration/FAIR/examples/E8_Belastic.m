%==============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: discretized elastic energy operator
%
%  - builds matrix B = getElasticMatrixStg
%  - compares to matrix-free implementation 
% see also getElasticMatrixStg and mfBy
%==============================================================================

clear, close all, help(mfilename);

% set parameter for the matrix based version
omega  = [0,1,0,2,0,3];
m      = [4,5,6];
mu     = 1;  
lambda = 1;
alpha  = 1;
B = getElasticMatrixStg(omega,m,mu,lambda);
yc = randn(size(B,2),1);
zc = randn(size(B,1),1);

figure(1); clf; spy(B); set(gca,'fontsize',30)
xlabel(['nz = ' int2str(nnz(B))],'fontsize',30);

% setup regularizer for the matrix free version
regularizer('reset','regularizer','mfElastic','alpha',alpha,'mu',mu,'lambda',lambda);
[Sc, dS, d2S] = regularizer(zeros(size(yc)),omega,m);

compare = @(s,l,r) fprintf('%-25s = %s\n',s,num2str(norm(l-r)));
compare('||B*yc-mfBy(y)||',B*yc,d2S.By(yc,omega,m));
compare('||B''*z-mfBy(z,''BTy'')||',B'*zc,d2S.BTy(zc,omega,m));

fprintf('<%s> done!\n',mfilename); pause;


