ang=[1:360];
N=50;
numThetan=30;
W_test=zeros(N,N,length(ang)*N,2);
sino=zeros(N,numThetan,length(ang)*N,2);
for slice = 1:N
    for an=1:length(ang)
        angle=ang(an);
        trans=1;
        do_setup;
        sino(:,:,an+length(ang)*(slice-1),1)=DisR_Simulated;
        W_test(:,:,an+length(ang)*(slice-1),1)=W;
        trans=0;
        do_setup;
        sino(:,:,an+length(ang)*(slice-1),2)=DisR_Simulated;
        W_test(:,:,an+length(ang)*(slice-1),2)=W;

    end
end
sino=reshape(sino,numThetan*(nTau+1),length(ang)*N*2);
W=reshape(W_test,N*N,length(ang)*N*2);
save('ML_train_36000.mat','sino','W');

