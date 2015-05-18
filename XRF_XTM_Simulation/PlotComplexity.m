A1=[...
25       1       0.193635      2.316696
100       1       0.483618      3.352679
400       1       7.711324      20.650408
1600       1       43.911890      295.759373
200       2       1.739629      35.223946
800       2       16.197967      563.913147
1600       4       37.332138      1129.284628
3200       2       119.428006      8839.254406];
A2=[...
25       1       0.198030      2.304490
100      1       1.272278      9.250259
400      1       16.028315      38.162054
1600       1       86.975281      521.684550
200       2       2.360553      69.492924
800       2       33.314572      1086.257282
1600       4       63.729197      2233.780828
3200       2       214.165188      17436.166420];
close all
n=4;
% semilogy(A1(1:n,1),A1(1:n,3),'rs-',A1(1:n,1),A1(1:n,4),'bs-',A2(1:n,1),A2(1:n,4),'m*-',A2(1:n,1),A2(1:n,3),'g*-', ...
% A1(n+1:end,1),A1(n+1:end,3),'cs--',A1(n+1:end,1),A1(n+1:end,4),'bs--',A2(n+1:end,1),A2(n+1:end,4),'m*--',A2(n+1:end,1),A2(n+1:end,3),'g*--' ...
% ,'linewidth',1.5);
% 
% legend('Analytic, 1 Angle, 1 Element', 'AD, 1 Angle (reverse), 1 Element', 'AD, 2 Angles (reverse), 1 Element','Analytic, 2 Angles, 1 Element' ...
% ,'Analytic, 1 Angle, more than 1 Element', 'AD, 1 Angle (forward), more than 1 Element', 'AD, 2 Angles (forward), more than 1 Element','Analytic, 2 Angles, more than 1 Element' ...
% ,'FontSize',12,'Location','NorthWest')

semilogy(A2(1:n,1),A2(1:n,4),'m*-',A2(1:n,1),A2(1:n,3),'g*-', ...
A2(n+1:end,1),A2(n+1:end,4),'m*--',A2(n+1:end,1),A2(n+1:end,3),'g*--' ...
,'linewidth',1.5);

hleg=legend( 'AD, 2 Angles (reverse), 1 Element','Analytic, 2 Angles, 1 Element' ...
, 'AD, 2 Angles (forward), more than 1 Element','Analytic, 2 Angles, more than 1 Element' ...
,'FontSize',12,'Location','Best');
set(hleg,'FontSize',16,'FontWeight','bold','Location','Best');
xlabel('$|\mathcal{V}||\mathcal{E}|$','FontSize',16,'FontWeight','bold','Interpreter','latex')
ylabel('Time Elapsed (sec)','FontSize',14,'FontWeight','bold')
% close all
% load ComA
% A=ComA;
% m=unique(A(:,1));
% rng('default');
% cmap = colormap(rand(length(m),3));
% for i=1:length(m); 
% ind=find(A(:,1)==m(i));
% semilogy(A(ind,2),A(ind,3),'.-','linewidth',2,'color',cmap(i,:));
% xlabel('Number of Angles','FontSize',16,'FontWeight','bold')
% ylabel('Time Elapsed (Log)','FontSize',14,'FontWeight','bold')
% axis([0 12 0 max(A(:,3))+100]);
% text(A(ind(end),2)+0.2,A(ind(end),3),['N=',num2str(m(i))],'FontSize',16,'FontWeight','bold')
% hold on;
% end

