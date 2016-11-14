Mt=-log(DisR_Simulated./I0);
Mtemp=Mt;
Mtemp(:,31:end)=Mtemp(end:-1:1,31:end);
Mt1=Mtemp;
for i=0:29,[cor,lag]=xcorr(Mtemp(:,1+i),Mtemp(:,31+i));[~,dd]=max(cor);d(i+1)=lag(dd);
    if(d(i+1)<=0)
        shift=abs(d(i+1))/2;
        Mt1(1:end-shift,i+1)=Mtemp(shift+1:end,i+1);
        Mt1(end-shift:end,i+1)=Mtemp(1:shift+1,i+1);
        Mt1(1:shift+1,i+31)=Mtemp(end-shift:end,i+31);
        Mt1(shift+1:end,i+31)=Mtemp(1:end-shift,i+31);
    else
        shift=d(i+1)/2;
        Mt1(1:end-shift,i+31)=Mtemp(shift+1:end,i+31);
        Mt1(end-shift:end,i+31)=Mtemp(1:shift+1,i+31);
        Mt1(1:shift+1,i+1)=Mtemp(end-shift:end,i+1);
        Mt1(shift+1:end,i+1)=Mtemp(1:end-shift,i+1);
    end
end

Mt1(:,31:end)=Mt1(end:-1:1,31:end); 
