for slice=2:201
    opt;
end
% beta=[1 0.0; 1 0.5; 1 1; 1 2; 1 6; 1 12; 1 24; 1 48; 1 100; 1 1e3];
% beta_dT=[1e2];
% beta=[1 0.5; 1 1;  1 6; 1 12; 1 24; 1 48; 1 100];
% for bt2=1:length(beta_dT)
%     beta_d=beta_dT(bt2);
% f_pareto=zeros(size(beta,1),2);
% for bt=1:size(beta,1)
%     TempBeta=beta(bt,1); Beta=beta(bt,2);
%     if(bt==1) 
%         initialize=1; 
%     else 
%         initialize=0;
%     end
%     opt;
%     [~,~,f_xrf,f_xrt]=feval(fctn,xstar);
%     f_pareto(bt,:)=[f_xrf,f_xrt];
% end
% save(['f_pareto',num2str(beta_d),'.mat'],'f_pareto');

