Nlevel=[3 6 12 24 48 96];
fid = fopen('Time_full_admm.txt','a');
for i=1:length(Nlevel)
    N=Nlevel(i);
    do_setup;
    optXRF;
    fprintf(fid,'%d       %f       %f\n',N^2*NumElement,t1,t2);
end
fclose(fid)

