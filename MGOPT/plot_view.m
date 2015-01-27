function plot_view(v,e)
%--------------------------------------------
% Plot the current search direction, and
% plot the step from the current v to v*
%--------------------------------------------
global vstar
global X N
global current_n

n = length(v);
j = find(N==current_n);

x     = X(1:N(j),j);
xstar = X(1:N(1),1);
vs    = interp1(xstar,vstar,x);

figure(12)

tt1 = sprintf('vs-v: Level = %2i n = %2i',j,current_n);
tt2 = sprintf('step: Level = %2i n = %2i',j,current_n);

subplot(1,3,1)
  plot(x,vs-v,'x',x,vs-v)
  title(tt1)
%  axis([0 1 -.2 .8])
subplot(1,3,2)
  plot(x,e,'x',x,e)
  title(tt2)
%  axis([0 1 -.2 .8])
subplot(1,3,3)
  r = (vs-v)./e;
  plot(x,r,'x',x,r)
  title('Ratio of Vectors')

disp('Hit any key to continue')
pause