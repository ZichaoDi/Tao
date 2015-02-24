% [M,I0,Ltol,thetan,omega,m,dTau]=InterLength;
% save('data/CornerSphere10theta90.mat','M','I0','Ltol','thetan','omega','m','dTau');
% return;
more off;
load('data/DisCornerSphere20theta90.mat');
fctn=@(I)sfun_discrete(I,M,I0,Ltol,thetan,omega,m,dTau);
x0=I0(:)+1*rand(size(I0(:)));%ones(prod(size(I0)),1);%
options = optimset('GradObj','on','Display','iter');
% [x,fval] = fminunc(fctn,x0,options);
[x,f,g] = tn(x0,fctn);
err=norm(x-I0(:));
fprintf('residule is %d \n',err);
figure('name','minimum')
imagesc(reshape(x,size(I0))); axis xy

% % report 10x10:
% %                                 Norm of      First-order 
% %  Iteration        f(x)          step          optimality   CG-iterations
% %      0        2.12211e+07                      2.58e+04                
% %      1        1.89001e+07             10       2.45e+04           2
% %      2        1.47026e+07             20       2.18e+04           2
% %      3        8.05704e+06             40       1.63e+04           2
% %      4        1.20021e+06             80       6.89e+03           2
% %      5            20398.4        69.6957            718           2
% %      6            226.966        14.8306           56.4           4
% %      7            4.23395        3.15128           6.61           7
% %      8          0.0238309       0.557445          0.588           9
% %      9        0.000210938      0.0353798         0.0491           9
% %     10        4.38741e-07     0.00362172        0.00334          10
% %     11        6.37898e-09    0.000118139       0.000272           7
% % 
% % Local minimum possible.
% % 
% % fminunc stopped because the final change in function value relative to 
% % its initial value is less than the default value of the function tolerance.
% % 
% % <stopping criteria details>
% % 
% % residule is 2.219341e-05 
%%%%%%====================================================================
% % 50x50
% %                                Norm of      First-order 
% % Iteration        f(x)          step          optimality   CG-iterations
% %     0        1.95128e+12                      4.93e+06                
% %     1         1.9495e+12             10       4.93e+06           2
% %     2        1.94594e+12             20       4.92e+06           2
% %     3        1.93882e+12             40       4.91e+06           2
% %     4        1.92462e+12             80        4.9e+06           2
% %     5         1.8964e+12            160       4.86e+06           2
% %     6        1.84064e+12            320       4.79e+06           2
% %     7        1.73185e+12            640       4.66e+06           2
% %     8         1.5252e+12           1280       4.38e+06           2
% %     9        1.15524e+12           2560        3.9e+06           2
% %    10        5.87364e+11           5100        2.9e+06           2
% %    11        2.33657e+11           5100       1.98e+06           2
% %    12        6.20877e+10           5100        1.1e+06           2
% %    13        5.34705e+09           5100       2.88e+05           3
% %    14        1.15612e+08        2469.98       3.57e+04           3
% %    15        2.31526e+06        908.828        2.6e+03           7
% %    16            27422.3        332.741            164          20
% %    17            463.815        68.0534           18.2          27
% %    18            1.64647        11.6226           1.15          35
% %    19          0.0152751       0.580804          0.098          35
% %    20        5.17553e-05      0.0627866        0.00664          39
% %    21         5.9436e-07     0.00260729       0.000601          30
% %    22        3.45778e-09    0.000432789       6.39e-05          36
% % 
% % Local minimum possible.
% % 
% % fminunc stopped because the final change in function value relative to 
% % its initial value is less than the default value of the function tolerance.
% % 
% % <stopping criteria details>
% % 
% % residule is 2.232288e-05 