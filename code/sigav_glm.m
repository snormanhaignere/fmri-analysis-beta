function MAT_file = sigav_glm(data_matrix_file, para_file, parameter_file, ...
    MAT_file, varargin)

% A GLM analysis is applied to signal averaged responses.
%
% Calculates percent signal change for a set of events/stimuli by signal
% averaging a fixed number of time-points after event onset. A set of predictors
% is then regressed against the signal averaged values, and the beta weights
% from this analysis are contrasted. Stats are computed using vanilla OLS
% equations, assuming independent, Gaussian errors, and a permutation test is
% used to compute stats by shuffling the order of conditions.
%
% 2016-09-09: Last edited, Sam NH
%
% 2016-09-21: Modified to deal with zero regressors and contrasts, Sam NH
%
% 2016-12-21: Made it possible to whiten the data before performing regression
% analysis. Doing so causes the regression analysis to perform LCMV.
%
% 2017-01-12: Moved some functionality to other helper functions that are shared
% by sigav_data.
%
% 2017-02-04: Saves z statistic for permutation test in addition to logP stat,
% Sam NH
%
% 2017-02-24: remove_unspecified_conditions no longer needed and removed
% 
% 2018-05-30: Add an optional tsnr threshold

% optional arguments and defaults
I.onset_delay = 5;
I.offset_delay = 1;
I.n_perms = 0;
I.whiten = false;
I.tsnr_threshold = 0;
I = parse_optInputs_keyvalue(varargin, I);

% analysis parameters
P = load(parameter_file);
P = glm_default_parameters(P);

%% Compute signal averaged responses

[psc, mean_signal, T, voxels_without_NaN, null_response] = sigav_data(...
    data_matrix_file, para_file, parameter_file, ...
    'onset_delay', I.onset_delay, 'offset_delay', I.offset_delay, ...
    'tsnr_threshold', I.tsnr_threshold);
n_trials = size(psc,1);
assert(length(T.onsets) == n_trials);

%% Regression and contrast matrix

[W, W_one_per_condition, C, ...
    nonzero_regressors, nonzero_contrasts, nonzero_conditions] = ...
    sigav_weights_and_contrasts(P, T);

%% Optionally whiten data and apply whitening matrix to the regressors as well

if I.whiten
    Z = (psc * psc')^(-1/2);
    psc = Z * psc;
    W = Z * W;
    W_one_per_condition = Z * W_one_per_condition;
end

%% Regression analyses

try
    % regress weighted event matrix
    % -> contrasts x voxels
    [beta_contrast, logP_ols, contrast_variance, df, residual] = ...
        regress_stats_ols(psc, W, C, 'demean', false); %#ok<ASGLU>
    
    % separate beta weight for regressor
    % redundant with previous analysis if contrast matrix is the identity
    % -> regressors x voxels
    beta_one_per_regressor = ...
        regress_stats_ols(psc, W, eye(size(W,2)), 'demean', false);
    
    % separate beta weight for each condition
    % -> condition x voxels
    beta_one_per_condition = ...
        regress_stats_ols(psc, W_one_per_condition, ...
        eye(size(W_one_per_condition,2)), 'demean', false);
    
catch error_message
    print_error_message(error_message);
    keyboard
    
end

%% Fill in NaN entries and save

% set NaN for zero contrasts
beta_contrast = fillin_NaN(beta_contrast, nonzero_contrasts, 1);
logP_ols = fillin_NaN(logP_ols, nonzero_contrasts, 1);
contrast_variance = fillin_NaN(contrast_variance, nonzero_contrasts, 1);
beta_one_per_regressor = fillin_NaN(beta_one_per_regressor, nonzero_regressors, 1);
beta_one_per_condition = fillin_NaN(beta_one_per_condition, nonzero_conditions, 1);

% set NaN for voxels with NaN
beta_contrast = fillin_NaN(beta_contrast, voxels_without_NaN, 2);
logP_ols = fillin_NaN(logP_ols, voxels_without_NaN, 2); %#ok<NASGU>
contrast_variance = fillin_NaN(contrast_variance, voxels_without_NaN, 2); %#ok<NASGU>
beta_one_per_regressor = fillin_NaN(beta_one_per_regressor, voxels_without_NaN, 2); %#ok<NASGU>
beta_one_per_condition = fillin_NaN(beta_one_per_condition, voxels_without_NaN, 2); %#ok<NASGU>
residual = fillin_NaN(residual, voxels_without_NaN, 2);
psc = fillin_NaN(psc, voxels_without_NaN, 2);
null_response = fillin_NaN(null_response, voxels_without_NaN, 2); %#ok<NASGU>
mean_signal = fillin_NaN(mean_signal, voxels_without_NaN, 2); %#ok<NASGU>

% save
save(MAT_file, 'psc', 'null_response', ...
    'mean_signal', 'voxels_without_NaN', 'beta_contrast', ...
    'logP_ols', 'contrast_variance', 'beta_one_per_condition', ...
    'beta_one_per_regressor', 'df', 'P', 'residual', '-v7.3');

%% Permutation test

% number of permuted weights
% optionally perform permutation test
if I.n_perms > 0
    
    [parent_directory,perm_MAT_file] = fileparts(MAT_file);
    perm_MAT_file = [parent_directory '/' perm_MAT_file '_' num2str(I.n_perms) 'perms.mat'];
    clear parent_directory;
    
    % estimate contrasts and residual from permutations
    beta_contrast_permtest = nan([I.n_perms, sum(nonzero_contrasts), sum(voxels_without_NaN)]);
    residual_permtest = nan([I.n_perms, sum(voxels_without_NaN)]);
    for i = 1:I.n_perms
        [beta_contrast_permtest(i,:,:), ~, ~, ~, residual_permtest(i,:)] = ...
            regress_stats_ols(psc(:,voxels_without_NaN), ...
            W(randperm(n_trials),:), C);
    end
    
    % convert to signed logP value
    [logP_permtest, zstat_permtest] = sig_via_null_gaussfit(...
        beta_contrast(nonzero_contrasts,voxels_without_NaN), ...
        beta_contrast_permtest);
    [logP_residual_permtest, zstat_residual_permtest] = sig_via_null_gaussfit(...
        residual(:,voxels_without_NaN), residual_permtest, 'tail', 'left');
    
    % set NaN for zero contrasts
    logP_permtest = fillin_NaN(logP_permtest, nonzero_contrasts, 1);
    zstat_permtest = fillin_NaN(zstat_permtest, nonzero_contrasts, 1);
    beta_contrast_permtest = fillin_NaN(beta_contrast_permtest, nonzero_contrasts, 2);
    
    % set NaN for voxels with NaN
    try
        logP_permtest = fillin_NaN(logP_permtest, voxels_without_NaN, 2); %#ok<NASGU>
        zstat_permtest = fillin_NaN(zstat_permtest, voxels_without_NaN, 2); %#ok<NASGU>
        beta_contrast_permtest = fillin_NaN(beta_contrast_permtest, voxels_without_NaN, 3); %#ok<NASGU>
        logP_residual_permtest = fillin_NaN(logP_residual_permtest, voxels_without_NaN, 2); %#ok<NASGU>
        zstat_residual_permtest = fillin_NaN(zstat_residual_permtest, voxels_without_NaN, 2); %#ok<NASGU>
        residual_permtest = fillin_NaN(residual_permtest, voxels_without_NaN, 2); %#ok<NASGU>
    catch me
        print_error_message(me);
        keyboard;
    end
    
    % save results to matfile
    save(perm_MAT_file, 'beta_contrast_permtest', 'logP_permtest',...
        'residual_permtest', 'logP_residual_permtest',...
        'zstat_permtest', 'zstat_residual_permtest');
    
end
