shift=xCOR(1:N_delta);
[~,~,aligned1]=feval(fctn_COR,xCOR);
for i=1:numThetan
    if(shift(i)>=0)
        aligned1(i,1:shift(i))=0;
    else
        aligned1(i,end+shift(i)-1:end)=0;
    end
end

