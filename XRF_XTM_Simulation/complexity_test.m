global N numThetan;
TimeS=[];

N=35;

%while i_t<=150
for i=1:10
    numThetan=i;
    m=[N,N];
    optXTM_XRF;
    TimeS=[TimeS;[prod(m),i,elapsedTime ]];
save TimeS TimeS
end
%    i_t=i_t+20;
%end
