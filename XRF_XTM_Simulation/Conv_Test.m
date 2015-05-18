%% Convergence Factor Test
%% Phantom Test (3 chemical elements) with 1 angle
do_setup;
fid = fopen('Convergence_Factor.txt','a');
for i=1:length(N)
    current_n=N(i);
    Joint=1;
    optXTM_XRF_Tensor;
    Joint=0;
    optXTM_XRF_Tensor;
    ntol=current_n^2*NumElement;
    fprintf(fid,'%d    %12.4e     %12.4e    %12.4e     %12.4e    %12.4e     %12.4e\n', ntol, convFac_J, t_J, errTol_J, convFac_XRF, t_XRF, errTol_XRF);
end

fclose(fid);

