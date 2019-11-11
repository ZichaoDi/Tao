% TV based image restoration using TwIST
% cameraman, blur uniform 9*9, BSNR = 40 dB
%
% This demo illustrates the computation of  
%
%     xe = arg min 0.5 ||Ax-y||^2 + tau TV(x)
%             x
% 
% with the TwIST algorithm. 

% set BSNR
more off;
BSNR = 40;
y=Mt0;
y=reshape(y,numThetan*(nTau+1),nchannel);
Py = var(y(:));
sigma= sqrt((Py/10^(BSNR/10)));

% smoothing parameter (empirical setting)
tau = 2e-2*sigma^2/0.56^2;

% extreme eigenvalues (TwIST parameter)
lam1=1e-4;    
% TwIST is not very sensitive to this parameter
% rule of thumb: lam1=1e-4 for severyly ill-conditioned% problems
%              : lam1=1e-1 for mildly  ill-conditioned% problems
%              : lam1=1    when A = Unitary matrix
% handle functions for TwIST
A=L;
tv_iters = 5;
% TV regularizer;
Phi = @(x) TVnorm(x);
varx = var(y(:)); 	% approximate var of x
x0 = x_base;%real(fft2(KC./(abs(KC).^2+10*sigma^2/varx).*ifft2(y)));
tolA = 1e-4;
% -- TwIST ---------------------------
% stop criterium:  the relative change in the objective function 
% falls below 'ToleranceA'
[x_twist,dummy,obj_twist,...
    times_twist,dummy,mse_twist]= ...
         TwIST(y,A,tau,...
         'lambda',lam1,...
         'Phi',Phi, ...
         'Monotone',1,...
         'Initialization',x0,...
         'StopCriterion',1,...
       	 'ToleranceA',tolA,...
         'MaxiterA',maxiter,...
         'Verbose', 1);

% figure(3); colormap gray; 
% title('TwIST restored image')
% imagesc(reshape(x_twist,N,N)); axis off;
% drawnow;

return;

     
% -- IST (lam1=1) ---------------------------
% stop criterium:  the objective function 
% falls below obj_twist(end)
%
% IST takes too long and thus we run only 200 iterations
[x_ist,dummy,obj_ist,...
  times_ist,dummy,mse_ist]= ...
         TwIST(y,A,tau,...
         'AT', AT, ...
         'lambda',1,...
         'True_x', x,...       
         'Psi', Psi, ...
         'Phi',Phi, ...
         'Monotone',1,...
         'MaxiterA',200,...
         'Initialization',x0,...
         'StopCriterion',3,...
       	 'ToleranceA',obj_twist(end),...
         'Verbose', 1);
     
     
     
figure(4); colormap gray; 
title('IST restored image')
imagesc(x_ist); axis off;
drawnow;
     
   
figure(5)
subplot(2,1,1)
semilogy(times_twist,obj_twist,'r',...
         times_ist,obj_ist,'b','LineWidth',2)
legend('TwIST','IST')
st=sprintf('tau = %2.2e, sigma = %2.2e',tau,sigma),...
title(st)
ylabel('Obj. function')
xlabel('CPU time (sec)')

grid
subplot(2,1,2)
plot(times_twist(2:end),mse_twist(2:end),'r',...
         times_ist(2:end),mse_ist(2:end),'b','LineWidth',2)
legend('TwIST','IST')
ylabel('MSE')
xlabel('CPU time (sec)')


fprintf(1,'TwIST   CPU time - %f\n', times_twist(end));
fprintf(1,'IST     CPU time - %f\n', times_ist(end));
