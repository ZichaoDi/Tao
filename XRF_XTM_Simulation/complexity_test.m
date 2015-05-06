global m;
% TimeS=zeros(50,2);
i_t=76;
while i_t<=150
    m=[i_t+1,i_t+1];
    optXRF;
    TimeS=[TimeS;[prod(m)*4,t]];
    i_t=i_t+20;
end