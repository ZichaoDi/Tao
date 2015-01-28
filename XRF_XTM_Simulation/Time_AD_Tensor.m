%% compare computation complexity between AdiMat and Tensor (Wendy's)
global current_n NumElement numThetan
TestSize=[...
          3 3 1 
          20 20 2
          40 40 2
          20 20 4
          40 40 4];
%           40 40 8];
fid = fopen('Time_AD_Tensor_1_2.txt','a'); 
fprintf(fid,'m_x    m_y      Ne      Tensor      AdiMat(secs)\n',numThetan);
for it=1:2
    numThetan=it;
    fprintf(fid,'Number of angles is %d =============\n',numThetan);
for i_size=1:size(TestSize,1)
    current_n=TestSize(i_size,1);
    NumElement=TestSize(i_size,3);
    optXTM_XRF_Tensor;
    fprintf(fid,'%d       %d       %d       %f      %f\n',TestSize(i_size,:),T_Tensor, T_AD); 
end 
end

fclose(fid);
    
