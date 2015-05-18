  function [xnew,Fnew,wnew,Jnew] = GNStep_VL(FName,JName,xc,Fc,wc,Jc,plist,Lmax,auto)
% [xnew,Fnew,wnew,Jnew] = GNStep(FName,JName,xc,Fc,wc,Jc,plist,Lmax,auto)
% Generates a Gauss-Newton Step
% 
%  FName     the name of a function F(x,plist) that accepts  a column
%            n-vectors and column m-vector
%  JName     the name of a function JF(x,plist) that returns the
%            Jacobian of F at x.
%     xc     a column n-vector, an approximate minimizer.         
%     Fc     F(xc,plist)
%     wc     Fc'*Fc/2
%     Jc     JF(xc,plist)
%  plist     parameter list
%   auto     0 for manual line search and nonzero otherwise.
%
%   xnew     a column n-vector, an improved minimizer.
%   Fnew     F(xnew,plist)
%   wnew     Fnew'*Fnew/2
%   Jnew     JF(xnew,plist)

   nplotvals = 20; % Number of F evaluations per line search.

% Get the Gauss-Newton  Search Direction

sc = -Jc\Fc;

% Line Search

% Try to get L<=Lmax so xc+L*sc is at least as good as xc.

L = Lmax;
Halvings = 0;
FL = feval(FName,xc+L*sc,plist);
while ((FL'*FL/2)>=wc) & Halvings<=10
   L = L/2;
   Halvings = Halvings+1;
   FL = feval(FName,xc+L*sc,plist);
end  
lambdavals = linspace(0,L,nplotvals);

wvals = zeros(nplotvals,1);
for k=1:nplotvals
   F = feval(FName,xc+lambdavals(k)*sc,plist);
   wvals(k) = (F'*F)/2;
end
if auto==0
   % Manual line search
   plot(lambdavals,wvals);
   xlabel('lambda');
   ylabel('w(xc+lambda*sc)');
   title('Enter the Best lambda value.');
   [lambda,y] = ginput(1);
   xnew = xc+lambda*sc;  
else
   % Automated line search
   
   [fnew,i] = min(wvals(1:nplotvals));
   xnew = xc + lambdavals(i(1))*sc;
end
Fnew = feval(FName,xnew,plist); 
wnew = .5*(Fnew'*Fnew);
Jnew = feval(JName,xnew,plist);