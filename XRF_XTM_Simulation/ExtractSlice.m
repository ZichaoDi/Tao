MAPS=[];
PWD=[pwd,'/data/ApsDataExtract/hong2_SV/'];
channels=dir([PWD,'/aligned_projections/']);
NumThetan=73;
channel_names=[];
for slice=1:51
    for i=3:length(channels)
        sub_file=dir(fullfile(PWD,'/aligned_projections/',channels(i).name));
        for thetan=1:NumThetan
            A=imread(fullfile(PWD,'/aligned_projections/',channels(i).name,sub_file(thetan+2).name));
            MAPS(i-2,thetan,:)=reshape(A(slice,:),1,1,size(A,2));
        end
        channel_names{i-2}=channels(i).name;
    end
    save(fullfile(PWD,['hong',num2str(slice)]),'MAPS','channel_names');
end

% slice=30;
% nTau=1750;
% numChannel=2000;
% spectra=zeros(nTau,NumThetan,numChannel);
% for t=1:NumThetan
%     data_t=h5read([PWD,'aligned_mca/raw/',num2str(172+t),'.h5'],'/raw');
%     spectra(:,t,:)=data_t(:,slice,:);
% end
% save([PWD,'spectra',num2str(slice),'.mat'],'spectra');
