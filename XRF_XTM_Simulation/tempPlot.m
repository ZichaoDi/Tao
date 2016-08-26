do_setup_mg;
d=xtm_level{1};
L=L_level{1};
return;
tic;x=lsqnonneg(L,d);time=toc; save('temp.mat','x','time');
return;
hollow=1;
NoSelfAbsorption=0;
do_setup_simulate;
save('XRF_rod_hollow.mat','XRF','DisR');
hollow=0;
NoSelfAbsorption=0;
do_setup_simulate;
save('XRF_rod_solid.mat','XRF','DisR');
NoSelfAbsorption=1;
do_setup_simulate;
save('XRF_rod_solid_noself.mat','XRF','DisR');
return;
global numThetan N Beta
Beta=1;
fid = fopen('Complexity_Tensor.txt','a');
for ii=15:5:30
    N=ii;
fprintf(fid,'=========================== %d\n', N);
for i=1:10
    numThetan=i;
    do_setup;
    optXTM_XRF_Tensor;
    fprintf(fid,'%d    %d      %f\n',prod(m)*NumElement,numThetan,T);
end
end
fclose(fid);
    



% close all
% figure(25);
% %  load test10
% 
%     clims=[0 max([xinitial;ws;xstar])];
%     for i=1:NumElement
%         subplot(4,NumElement,i);
%         
%         errCom=reshape(xinitial(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-ws(prod(m)*i-prod(m)+1:prod(m)*i
%         imagesc(errCom,clims);colormap gray
%         if(i==1)
%             ylabel('Initial Guess','fontsize',12)
%         end
%         title(Element{Z(i)},'fontsize',12);
%     end
%     
%     for i=1:NumElement
%         subplot(4,NumElement,i+1*NumElement);
%         
%         errCom=reshape(xstar(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));%-ws(prod(m)*i-prod(m)+1:prod(m)*i
%         imagesc(errCom,clims);colormap gray
%         if(i==1)
%             ylabel('Final Soluction','fontsize',12)
%         end
%     end
%     
%     for i=1:NumElement
%         subplot(4,NumElement,i+2*NumElement);
%         
%         errCom=reshape(ws(prod(m)*i-prod(m)+1:prod(m)*i),m(1),m(2));
%         imagesc(errCom,clims);colormap gray
%         if(i==1)
%             ylabel('True Soluction','fontsize',12)
%         end
%         if(i==NumElement)
%         hp4 = get(subplot(4,NumElement,i+2*NumElement),'Position');
%         colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(3)]);
%         end
%     end
%     for i=1:NumElement; subplot(4,NumElement,i+3*NumElement);
%         plot(1:prod(m),sort(xinitial(prod(m)*i-prod(m)+1:prod(m)*i)),'ro',1:prod(m),sort(xstar(prod(m)*i-prod(m)+1:prod(m)*i)),'bs',1:prod(m),sort(ws(prod(m)*i-prod(m)+1:prod(m)*i)),'g*')
%         xlim([0 prod(m)]);
%         if(i==1)
%             hleg=legend('initial','final','optimal','FontSize',3);
% %             set(hleg,'units','pixels');
% %             lp=get(hleg,'outerposition');
% %             set(hleg,'Location','NorthWest', 'Box', 'off');%,'outerposition',[lp(1)-30,80,lp(3),30]);
%             ylabel('solution','fontsize',12)
%         end
%     end
