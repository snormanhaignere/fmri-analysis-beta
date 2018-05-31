function MAT_file = sigav_inverse_glm(data_matrix_file, para_file, parameter_file, ...
    MAT_file, varargin)

% Signal averaged voxel responses are regressed against a desired response to
% learn a weighted average of the voxel responses that best approximates that
% response. Ridge is used to regularize the solution. 
% 
% 2017-01-12: Started coding, Sam NH
% 
% 2018-05-30: Added optional tSNR threshold

% optional arguments and defaults
I.onset_delay = 5;
I.offset_delay = 1;
I.n_perms = 0;
I.folds = 5;
I.method = 'ridge';
I.K = 2.^(-100:100);
I.std_feats = false;
I.demean_feats = false;
I.rand_seed_for_trial_order = 1;
I.tsnr_threshold = 0;
I = parse_optInputs_keyvalue(varargin, I);

% analysis parameters
P = load(parameter_file);
P = glm_default_parameters(P);

%% Compute signal averaged responses

[psc, mean_signal, T, voxels_without_NaN] = sigav_data(...
    data_matrix_file,  para_file, parameter_file, ...
    'onset_delay', I.onset_delay, 'offset_delay', I.offset_delay, ...
    'tsnr_threshold', I.tsnr_threshold);
n_trials = size(psc,1);
assert(length(T.onsets) == n_trials);

%% Regressors and contrast matrix

[W, ~, C, nonzero_regressors, nonzero_contrasts] = sigav_weights_and_contrasts(P, T);

%% Regression analyses

n_betas = size(W,2);
n_voxels_without_NaN = sum(voxels_without_NaN);
beta_one_per_regressor = nan(n_betas, n_voxels_without_NaN);
offset = nan(n_betas,1);
for i = 1:size(W,2)
    try
        ResetRandStream2(I.rand_seed_for_trial_order);
        xi = randperm(n_trials);
        x = regress_weights_from_2way_crossval(...
            psc(xi,:), W(xi,i), 'folds', I.folds, ...
            'method', I.method, 'K', I.K, ...
            'std_feats', I.std_feats, ...
            'demean_feats', I.demean_feats);
        beta_one_per_regressor(i,:) = x(2:end);
        offset(i) = x(1);
        clear x;
    catch ME
        print_error_message(ME);
        keyboard;
    end
end
clear xi;

% contrasts
beta_contrast = C' * beta_one_per_regressor;

%% Fill in NaN entries and save

% set NaN for zero contrasts
beta_contrast = fillin_NaN(beta_contrast, nonzero_contrasts, 1);
beta_one_per_regressor = fillin_NaN(beta_one_per_regressor, nonzero_regressors, 1);
offset = fillin_NaN(offset, nonzero_regressors, 1); %#ok<NASGU>

% set NaN for voxels with NaN
beta_contrast = fillin_NaN(beta_contrast, voxels_without_NaN, 2); %#ok<NASGU>
beta_one_per_regressor = fillin_NaN(beta_one_per_regressor, voxels_without_NaN, 2); %#ok<NASGU>
psc = fillin_NaN(psc, voxels_without_NaN, 2); %#ok<NASGU>
mean_signal = fillin_NaN(mean_signal, voxels_without_NaN, 2); %#ok<NASGU>

% save
save(MAT_file, 'psc', 'mean_signal', 'voxels_without_NaN', 'beta_contrast', ...
    'beta_one_per_regressor', 'offset', 'P', '-v7.3');

%% Permutation test

% number of permuted weights
% optionally perform permutation test
if I.n_perms > 0
    
    error('Permutation test not setup yet.');
    
    %     [parent_directory,perm_MAT_file] = fileparts(MAT_file);
    %     perm_MAT_file = [parent_directory '/' perm_MAT_file '_' num2str(I.n_perms) 'perms.mat'];
    %     clear parent_directory;
    %
    %     % estimate contrasts and residual from permutations
    %     beta_contrast_permtest = nan([I.n_perms, sum(nonzero_contrasts), sum(voxels_without_NaN)]);
    %     residual_permtest = nan([I.n_perms, sum(voxels_without_NaN)]);
    %     for i = 1:I.n_perms
    %         [beta_contrast_permtest(i,:,:), ~, ~, ~, residual_permtest(i,:)] = ...
    %             regress_stats_ols(psc(:,voxels_without_NaN), ...
    %             W(randperm(n_trials),:), C);
    %     end
    %
    %     % convert to signed logP value
    %     logP_permtest = sig_via_null_gaussfit(...
    %         beta_contrast(nonzero_contrasts,voxels_without_NaN), ...
    %         beta_contrast_permtest);
    %     logP_residual_permtest = sig_via_null_gaussfit(...
    %         residual(:,voxels_without_NaN), residual_permtest, 'tail', 'left');
    %
    %     % set NaN for zero contrasts
    %     logP_permtest = fillin_NaN(logP_permtest, nonzero_contrasts, 1);
    %     beta_contrast_permtest = fillin_NaN(beta_contrast_permtest, nonzero_contrasts, 2);
    %
    %     % set NaN for voxels with NaN
    %     logP_permtest = fillin_NaN(logP_permtest, voxels_without_NaN, 2); %#ok<NASGU>
    %     beta_contrast_permtest = fillin_NaN(beta_contrast_permtest, voxels_without_NaN, 3); %#ok<NASGU>
    %     logP_residual_permtest = fillin_NaN(logP_residual_permtest, voxels_without_NaN, 2); %#ok<NASGU>
    %     residual_permtest = fillin_NaN(residual_permtest, voxels_without_NaN, 2); %#ok<NASGU>
    %
    %     % save results to matfile
    %     save(perm_MAT_file, 'beta_contrast_permtest', 'logP_permtest',...
    %         'residual_permtest', 'logP_residual_permtest');
    
end
