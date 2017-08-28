function [Gsim, similarity_metric_MAT_file] = ...
    noisecorr_similarity_via_regression_voxelgrid(...
    G, F, output_directory, varargin)

% Takes a feature matrix F (stim x feature) and voxel data structure G in the
% form of a grid, and computes a noise-corrected measure of the similarity
% of the response against the a prediction derived from cross-validated linear
% regression. Divides the analysis up into batches and makes it possible to use
% slurm to parralize the analysis on openmind. 
% 
% 2017-03-18: Added similarity_metric to the file name to which the similarity
% measures are saved.
%
% 2017-08-11: Making within-fold correlation optional

if ischar(F)
    load(F, 'F');
end

if ischar(G)
    load(G, 'G');
end

I.test_folds = 2;
I.train_folds = min(max(round((size(F,1) / I.test_folds) / 10),2),10);
I.similarity_metric = 'pearson';
I.regression_method = 'ridge';
I.K = [];
I.slurm = false;
I.overwrite_predictions = false;
I.overwrite_correlations = false;
I.mem = '8000';
I.std_feats = true;
I.demean_feats = true;
I.batch_size = 1000;
I.groups = [];
I.n_perms = 0;
I.correction_method = 'correlation-based';
I.regularization_metric = 'unnormalized-squared-error';
I.noisecorr_regmetric = false;
I.max_run_time_in_min = num2str(60*10);
I.within_fold_correlation = true;
I.only_cross_column_corr = false;
I = parse_optInputs_keyvalue(varargin, I);

% cannot use correlation-based correction for variance-based metric
if strcmp(I.similarity_metric, 'demeaned-squared-error')
    I.correction_method = 'variance-based';
end

%% Predictions

prediction_MAT_file = [output_directory '/predictions.mat'];
if ~exist(prediction_MAT_file, 'file') || I.overwrite_predictions
    
    % stimuli x rep x voxel
    D = grid2matrix(G);
    [n_stim, n_rep, n_vox] = size(D); %#ok<ASGLU>
    
    % remove zero or NaN voxels
    nonempty_voxels = ~(all(all(D < 1e-100,1),2) | all(all(isnan(D),1),2));
    D = D(:,:,nonempty_voxels);
    n_nonempty_vox = size(D,3);
    
    % stim x (rep * voxel)
    if ~I.noisecorr_regmetric
        D = reshape(D, n_stim, n_rep  * n_nonempty_vox);
    end
    
    % directory to save results for each voxel
    directory_to_save_batch_results = ...
        [output_directory '/batch_predictions'];
    if ~exist(directory_to_save_batch_results, 'dir')
        mkdir(directory_to_save_batch_results);
    end
    [Yh, folds] = regress_predictions_wrapper(...
        F, D, 'test_folds', I.test_folds, 'train_folds', I.train_folds, ...
        'method', I.regression_method, 'K', I.K, ...
        'output_directory', directory_to_save_batch_results, ...
        'mem', I.mem, 'save_results', true, 'slurm', I.slurm, ...
        'std_feats', I.std_feats, 'demean_feats', I.demean_feats, ...
        'groups', I.groups, ...
        'regularization_metric', I.regularization_metric, ...
        'noisecorr_regmetric', I.noisecorr_regmetric, ...
        'correction_method', I.correction_method, ...
        'batch_size', I.batch_size, ...
        'max_run_time_in_min', I.max_run_time_in_min, ...
        'overwrite', I.overwrite_predictions);
    
    save(prediction_MAT_file, 'Yh', 'D', 'folds', 'nonempty_voxels',...
        'n_stim', 'n_rep', 'n_vox', 'n_nonempty_vox');
    
end

%% Similarity metric

% noise-corrected correlations
similarity_metric_MAT_file = ...
    [output_directory '/noise_corrected_' I.similarity_metric '.mat'];

% delete file if is un-loadable
if exist(similarity_metric_MAT_file, 'file')
    try 
        load(similarity_metric_MAT_file, 'r');
    catch
        delete(similarity_metric_MAT_file);
    end
end

if ~exist(similarity_metric_MAT_file, 'file') || I.overwrite_correlations
    
    % load variables from previous step
    load(prediction_MAT_file, 'Yh', 'D', 'folds', 'nonempty_voxels', ...
        'n_stim', 'n_rep', 'n_vox', 'n_nonempty_vox');
    
    % separate out reps / voxels into different dimensions
    D = reshape(D, [n_stim, n_rep, n_nonempty_vox]);
    Yh = reshape(Yh, [n_stim, n_rep, n_nonempty_vox]);
    
    if ~I.within_fold_correlation
        folds = ones(size(folds));
    end
    
    % compute correlation
    r = nan(n_nonempty_vox,1);
    for i = 1:size(D,3);
        switch I.correction_method
            case 'correlation-based'
                r(i) = normalized_correlation_within_folds(...
                    D(:,:,i), Yh(:,:,i), folds, 'metric', I.similarity_metric, ...
                    'only_cross_column_corr', I.only_cross_column_corr);
            case 'variance-based'
                r(i) = noise_corrected_similarity_within_folds(...
                    D(:,:,i), Yh(:,:,i), folds, 'metric', I.similarity_metric,...
                    'only_cross_column_cov', I.only_cross_column_corr);
            otherwise
                error('Switch statement fell through');
        end
    end
    
    % add back in empty voxels as NaNs
    r = fillin_NaN(r, nonempty_voxels, 1);
    
    save(similarity_metric_MAT_file, 'r');
    
else
    
    load(similarity_metric_MAT_file, 'r');
    
end

%% Permutation test
% 
% if I.n_perms > 0
%     
%     % noise-corrected correlations
%     permtest_MAT_file = [output_directory '/permtest_' num2str(I.n_perms) 'perms' ...
%         '_' I.similarity_metric '.mat'];
%     if ~exist(permtest_MAT_file, 'file') || I.overwrite_correlations
%         
%         % load variables from previous step
%         load(prediction_MAT_file, 'Yh', 'D', 'folds', 'nonempty_voxels', ...
%             'n_stim', 'n_rep', 'n_vox', 'n_nonempty_vox');
%         
%         % separate out reps / voxels into different dimensions
%         D = reshape(D, [n_stim, n_rep, n_nonempty_vox]);
%         Yh = reshape(Yh, [n_stim, n_rep, n_nonempty_vox]);
%         
%         % compute denominator of similarity metric
%         sig = nan(n_nonempty_vox,1);
%         tic;
%         for i = [1771        1791        1814        1815        1935        2016        2018] %1:size(D,3);
%             
%             if mod(i,100) == 0
%                 toc;
%                 fprintf('Voxel %d\n', i);
%                 drawnow;
%                 tic;
%             end
%             
%             % denominator without permutation
%             [~,~,d] = normalized_correlation_within_folds_v2(...
%                 D(:,:,i), Yh(:,:,i), folds, 'metric', I.similarity_metric);
%             
%             % denominator with permutation
%             d_perms = nan(I.n_perms,1);
%             Dp = nan(size(D,1), size(D,2));
%             Yhp = nan(size(Yh,1), size(Yh,2));
%             for j = 1:I.n_perms
%                 for l = 1:size(D,2)
%                     perm_order = nan(n_stim,1);
%                     for k = 1:max(folds)
%                         xi = find(folds==k);
%                         perm_order(xi) = xi(randperm(length(xi)));
%                     end
%                     Dp(:,l) = D(perm_order,l,i);
%                 end
%                 
%                 for l = 1:size(Yh,2)
%                     perm_order = nan(n_stim,1);
%                     for k = 1:max(folds)
%                         xi = find(folds==k);
%                         perm_order(xi) = xi(randperm(length(xi)));
%                     end
%                     Yhp(:,l) = Yh(perm_order,l,i);
%                 end
%                 
%                 [~,~,d_perms(j)] = normalized_correlation_within_folds_v2(...
%                     Dp, Yhp, folds, 'metric', I.similarity_metric);
%             end
%             
%             % significance metric based on counts
%             if ~isnan(d)
%                 d_perms(isnan(d_perms)) = 0;
%                 sig(i) = sig_via_null_counts(d, d_perms);
%             end
%             
%         end
%         
%         % add back in empty voxels as NaNs
%         sig = fillin_NaN(sig, nonempty_voxels, 1);
%         
%         save(permtest_MAT_file, 'sig');
%         
%     else
%         
%         load(permtest_MAT_file, 'sig');
%         
%     end
%     
%     save(similarity_metric_MAT_file, '-append', 'sig');
% 
% end

%% Reformating

% convert similarity metric to grid
Gsim = G;
Gsim.grid_data{1} = nan(size(Gsim.grid_data{1},1), size(Gsim.grid_data{1},2));
Gsim.grid_data{2} = nan(size(Gsim.grid_data{2},1), size(Gsim.grid_data{2},2));
Gsim = matrix2grid(r, Gsim);
save(similarity_metric_MAT_file, '-append', 'Gsim');

% % convert significance metric to grid
% if I.n_perms > 0
%     Gsig = G;
%     Gsig.grid_data{1} = nan(size(Gsig.grid_data{1},1), size(Gsig.grid_data{1},2));
%     Gsig.grid_data{2} = nan(size(Gsig.grid_data{2},1), size(Gsig.grid_data{2},2));
%     Gsig = matrix2grid(sig, Gsig);
%     save(similarity_metric_MAT_file, '-append', 'Gsig');
% else
%     Gsig = [];
% end

