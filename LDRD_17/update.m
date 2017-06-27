function vh=update(vH,res_prob)
global N_level NumElement current_n
%== Full Weighting to interpolate
%== res_prob=0: upsample data; res_prob=1: upsample variable
j=find(N_level==current_n);
vH=reshape(vH,N_level(j),N_level(j),NumElement);
vh=restrict_residual(N_level(j-1))'*vH*restrict_residual(N_level(j-1));
vh=vh(:)/2;
