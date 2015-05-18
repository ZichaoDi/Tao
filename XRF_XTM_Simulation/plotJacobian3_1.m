% AX=[...
%     0 0    14 0 0  14  0 0 14
%     0 0   217 0 0  216 0 0 216
%     0 14    0 0 14  0  0 14 0
%     0 217 0 0 215 0 0 215 0
%     13 0 0 14 0 0 14 0 0
%     212 0 0 213 0 0 214 0 10]; %%Jacobian of XRF
% AT=[...
%     0 0 -31 -14 0 -31 0 0 -31
%     -13 -31 0 0 -31 0 0 -31 0
%     -31 -13 0 -31 -1 0 -31 0 0]; %%Jacobian of XRT

% load JacobImprove
% AX=JR;
% AT=JT;
p=zeros(9,4);
load J3_under;
AJ=[AX;1e8*AT];%%Jacobian of JRT



% GET SOME NULL SPACES and DIRECTIONS
NX = null(AX);
NT = null(AT);
DX = sum(NX,2);
DT = sum(NT,2);
p(:,1)=DT;
p(:,2)=DX;

load J3_over;
AJ=[AX;1e8*AT];%%Jacobian of JRT



% GET SOME NULL SPACES and DIRECTIONS
NX = null(AX);
NT = null(AT);
DX = sum(NX,2);
DT = sum(NT,2);

DX=DX/norm(DX)*1.5e-5;
DT=DT/norm(DT)*1.5e-5;

p(:,3)=DT;
p(:,4)=DX;


lim=1;%e-5;
close all
% a=zeros(9,1); 
% a(5)=-1;a(8)=1;
% p=a;



%%%===========================================================
% load Example3_1.mat;
% do_setup;

% close all
% 
% xs=ones(9,1);
%  load xsX3over_1;
%  xsX=xstar;
%  load xsT3over;
% xsT=xstar;
% % 
% % load xsJ3over;
% % xsJ=xstar;
% p=xsX-xs;
% % ind=[1 2 3 4 5 7 8];
% % p(ind)=0;
% lim=1;
% norm(AT*p)


t=[linspace(-lim,-eps,10),0,linspace(eps,lim,10)];

numRow=1;
for itt=4
    ft=[];
for i=1:length(t),
d=t(i)*p(:,itt);
%      f(i,1)=norm(AX*d);f(i,2)=norm(AT*d); f(i,3)=norm(AJ*d);
    [f]=feval(fctn,ones(9,1)+d);
    [f_xrf]=feval(fctn_X,ones(9,1)+d);
    f_xrt=feval(fctn_T,ones(9,1)+d);
    ft(i,1)=f_xrf;ft(i,2)=f_xrt; ft(i,3)=f;
end

plot(t,ft(:,1),'ro-',t,ft(:,2),'gs-',t,ft(:,3),'b.-')


% load Example3_1_2.mat
save('./result/MatFile/Example3_1_2R.mat','t','ft')
% plot(t,ft(:,1),'ro-',t,ft(:,2),'gs-',t,ft(:,3),'b.-')
axis([min(t) max(t) -1 max(f(:))])
FS = 16;
LW = 1.5;
MS = 10;
hr=legend('XRF','XRT','JRT');
set(findall(gcf,'type','line'),'LineWidth',LW,'MarkerSize',MS)

set(hr,'FontSize',FS,'FontWeight','bold','Location','northeast');%'southeast');
xlabel('Neighborhood of Local Minimum','FontSize',FS,'FontWeight','bold')
ylabel('Residual','FontSize',FS,'FontWeight','bold')

set(gca,'FontSize',FS-2)
axis tight
box on

print('-depsc', '~/Documents/TAO/TAOlocal/optdi/Reconstruction/Writing/figures/JacobianExpIllustrationOverT')
end

