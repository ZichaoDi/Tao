base_path = pwd;
increment = 1;
max_angle_range = 180;

%% load data and true shifts
load(fullfile(base_path,'projections360.mat'))
load(fullfile(base_path,'shifts_true360.mat'))
shifts_true = shifts;

output_path = base_path;
if ~exist(output_path,'dir'); mkdir(output_path); end

%% resample data
theta = -90:0.5:89.5; % rotation angles

switch max_angle_range
    case 180
        ind_start = 1;
        ind_end = 360;
    case 150
        ind_start = 31;
        ind_end = 330;
    case 120
        ind_start = 61;
        ind_end = 300;
end
theta = theta(ind_start:increment:ind_end);
shifts_true = shifts_true(:,ind_start:increment:ind_end,:);
projections = projections(:,:,ind_start:increment:ind_end);
Nangles = length(theta);

%% start alignment
N_test = 1; % number of tests
shifts_aligned = zeros(N_test,Nangles,2);

for k=1:N_test
    disp(k)

    %% Shift true projections to simulate misaligned data
    projections_unaligned = zeros(size(projections));
    for i=1:size(projections,3)
        projections_unaligned(:,:,i) = circshift(projections(:,:,i),[shifts_true(k,i,2),shifts_true(k,i,1)]);
    end

    %% use your favorite method to find out shifts using projections_unaligned


    %% EXAMPLE: Optical flow alignment
    zero_degree_index = find(theta==0);

    %% Find random affine transformation
    
    %optimizer = registration.optimizer.RegularStepGradientDescent;
    
    %metric = registration.metric.MattesMutualInformation;
    metric = registration.metric.MeanSquares;
    
    optimizer = registration.optimizer.OnePlusOneEvolutionary;
    optimizer.InitialRadius = 6.25e-3;
    optimizer.Epsilon = 1.5e-6;
    optimizer.GrowthFactor = 1.05;
    optimizer.MaximumIterations = 300;
    
    %optimizer.MaximumIterations = 1000;
    %optimizer.MinimumStepLength = 1e-6;
    %optimizer.GradientMagnitudeTolerance = 1e-5;
    
    for i=zero_degree_index+1:Nangles
        %disp(i)
        ref_img = projections_unaligned(:,:,i-1);
        mov_img = projections_unaligned(:,:,i);
        tform = imregtform(ref_img, mov_img, 'translation', optimizer, metric,'DisplayOptimization',false);
        shifts_temp = tform.T(3,1:2);
        shifts_aligned(k,i,1) = shifts_aligned(k,i-1,1) + shifts_temp(1);
        shifts_aligned(k,i,2) = shifts_aligned(k,i-1,2) + shifts_temp(2);
    end
    
    for i=zero_degree_index-1:-1:1
        %disp(i)
        ref_img = projections_unaligned(:,:,i+1);
        mov_img = projections_unaligned(:,:,i);
        tform = imregtform(ref_img, mov_img, 'translation', optimizer, metric,'DisplayOptimization',false);
        shifts_temp = tform.T(3,1:2);
        
        shifts_aligned(k,i,1) = shifts_aligned(k,i+1,1) + shifts_temp(1);
        shifts_aligned(k,i,2) = shifts_aligned(k,i+1,2) + shifts_temp(2);
    end
  
end
%% save retulsts
%{
save_name = sprintf('Np%d_increment%.1f_SNR_%d_OPOE_MSE.mat', size(projections,3), increment*0.5, SNR);
save_path = fullfile(base_path,'matlab_imregister/');
if ~isfolder(save_path); mkdir(save_path); end
save(strcat(save_path,save_name),'shifts_aligned');
%}

