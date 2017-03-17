function second_level_ols(MAT_files_first_level, MAT_file_second_level, varargin)

% 2016-08-26: Added residual statistics, Sam NH
%
% 2016-08-27: Modified how optional arguments are handled, current
% implementation requires the MAT files the results are saved to, to be
% specified, Sam NH
%
% 2016-09-09: Results of permutation test saved as a separate MAT file, Sam NH
% 
% 2016-09-09: No longer calls first level script, glm_event_regression.m, Sam NH
% 
% 2016-09-10: Made residual permutation test one-tailed, Sam NH
% 
% 2016-09-21: Modified to deal with NaN values, Sam NH
% 
% 2017-02-24: Removed data_matrix_files and para_files which weren't used
% 
% 2017-02-24: Removed permutation test, moved to second_level_permtest.m

% number of total runs
n_runs = length(MAT_files_first_level);
assert(n_runs == length(MAT_files_first_level));

%% Second level, ols stats

% load all runs
% beta_contrast_allruns: runs x contrast x voxel
for i = 1:n_runs
    
    % load first level analysis
    X = load(MAT_files_first_level{i},'beta_contrast',...
        'contrast_variance','P','df','residual','beta_one_per_regressor');
    
    if i == 1
        beta_contrast_allruns = nan([n_runs, size(X.beta_contrast)]);
        beta_one_per_regressor_allruns = nan([n_runs, size(X.beta_one_per_regressor)]);
        contrast_variance_allruns = nan([n_runs, size(X.contrast_variance)]);
        residual_allruns = nan([n_runs, size(X.residual)]);
        dfs = nan(n_runs,1);
        P = X.P;  %#ok<NASGU>
    end
    try
        beta_contrast_allruns(i,:,:) = X.beta_contrast;
        beta_one_per_regressor_allruns(i,:,:) = X.beta_one_per_regressor;
        contrast_variance_allruns(i,:,:) = X.contrast_variance;
        residual_allruns(i,:,:) = X.residual;
        dfs(i) = X.df;
    catch me
        print_error_message(me);
        keyboard
    end
end
clear X;

% set beta contrast and residual to the average across runs
beta_contrast = squeeze_dims( nanmean(beta_contrast_allruns, 1), 1); %#ok<NASGU>
beta_one_per_regressor = squeeze_dims( nanmean(beta_one_per_regressor_allruns, 1), 1); %#ok<NASGU>
residual = squeeze_dims( mean(residual_allruns, 1), 1); %#ok<NASGU>

% ols stats across runs
[logP_fixed, contrast_variance_fixed] = ...
    fixed_effects(beta_contrast_allruns, contrast_variance_allruns, dfs); %#ok<ASGLU>
[logP_random, contrast_variance_random] = ...
    random_effects(beta_contrast_allruns); %#ok<ASGLU>

% save results
save(MAT_file_second_level, 'beta_contrast', 'logP_fixed', 'logP_random', ...
    'contrast_variance_fixed', 'contrast_variance_random', 'P', 'residual', ...
    'beta_one_per_regressor');

% clear variables no longer needed
clear beta_contrast_allruns beta_one_per_regressor_allruns beta_one_per_regressor ...
    contrast_variance_allruns logP_fixed contrast_variance_fixed logP_random ...
    contrast_variance_random;
