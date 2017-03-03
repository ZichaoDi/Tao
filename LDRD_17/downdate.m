function vH=downdate(vh,res_prob)
global N_level current_n NumElement numThetan LevelScale
%== Full Weighting to restrict
%== res_prob => 0: downsample variable; 1: downsample gradient; 3: downsample data
j=find(N_level==current_n);
if(res_prob==3) %% vh: numThetan X nTau+1 X numChannel
    
    for n=1:numThetan
        if(ndims(vh)==2)
            vH(:,n)=restrict_residual(size(vh,1))*vh(:,n);
        elseif(ndims(vh)==3)
            vH(n,:,:)=restrict_residual(size(vh,2))*vh(n,:,:);
        end
    end
else
    vh=reshape(vh,N_level(j-1),N_level(j-1),NumElement);
    vH=restrict_residual(N_level(j-1))*vh*restrict_residual(N_level(j-1))';
    vH=vH(:)/1;
end
