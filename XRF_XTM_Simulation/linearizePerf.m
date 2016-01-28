fid = fopen('LinearizeDirection.txt','a');
fprintf(fid, '============================== pcg: # of variables is %d \n',num2str(N(1)));
fprintf(fid,'ss   sg    gg    angle(s,g)   f(1)    f(2)   f(3)    f(4)   f(0)\n');
alpha=1./(2.^(0:10));
alpha=[1 1e-1 1e-2 1e-3];
h_admm=[];
f1=[];
gg=[];
steep=[];
f_line=[];
for i=1:icycle;
    [f1(i),gg(:,i)]=feval(fctn0,x_admm(:,i));
    h_admm(:,i)=x_admm(:,i+1)-x_admm(:,i);
    for step=1:length(alpha)
        f_line(step,i)=feval(fctn0,x_admm(:,i)+alpha(step)*h_admm(:,i));
    end
    df(i)=h_admm(:,i)'*gg(:,i);
    alterD(i)=h_admm(:,i)'*h_admm(:,i);
    steep(i)=gg(:,i)'*gg(:,i);
    ang(i)=acosd(dot(h_admm(:,i),gg(:,i)));
    fprintf(fid,'%3.5e       %3.5e       %3.5e       %d      %e     %e      %e    %e     %e\n',alterD(i), df(i), steep(i), real(ang(i)), f_line(:,i), f1(i)); 
    i
end
[f1(icycle+1),gg(:,icycle+1)]=feval(fctn0,x_admm(:,end));
steep(icycle+1)=gg(:,end)'*gg(:,end);
res=f1(2:end)-f1(1:end-1);
figure, subplot(1,2,1),plot(1:icycle,df,'r.-',1:icycle,res,'bo-') 
subplot(1,2,2),semilogy(1:icycle,abs(df),'r.-',1:icycle,abs(res),'bo-',0:icycle,steep,'g*-',1:icycle,alterD,'ks-')
legend('|f-fold|','|transpose(x-xold)*g|','transpose(g)*g','transpose(x-xold)*(x-xold)')
xlabel('cycle #');
figure, for i=1:10; subplot(2,5,i);semilogy(alpha,f_line(:,i),'b*-');title(['cycle #',num2str(i)]);end
fclose(fid);

