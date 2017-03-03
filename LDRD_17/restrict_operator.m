function IhH=restric_operator(m); 
global N

j=find(N==m);
I_h_2h=zeros(N(j),N(j-1));
if(mod(N(j),2)~=0)
    for i =1:N(j)
        I_h_2h(i,2*i-1)=1/2;
        I_h_2h(i,2*i)=1;
        I_h_2h(i,2*i+1)=1/2;
    end    
    IhH=1/2*I_h_2h(:,2:end-1);
else
    for i =1:N(j)
        I_h_2h(i,2*i-1)=1;
        I_h_2h(i,2*i)=1;
    end   
    IhH=I_h_2h;
end
