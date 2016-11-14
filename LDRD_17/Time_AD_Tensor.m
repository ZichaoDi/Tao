%% compare computation complexity between AdiMat and Tensor (Wendy's)
global current_n NumElement numThetan
TestSize=[...
          5  5  1
          10 10 1 
          20 20 1
          40 40 1
          10 10 2
          20 20 2
          40 40 2
          20 20 4];
fid = fopen('Time_AD_Tensor_1_2.txt','a');
fprintf(fid,'If single element, use reverse mode of AD, otherwise, use forward mode\n');
fprintf(fid,'m_x    m_y      Ne      Tensor      AdiMat(secs)\n');
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
    
