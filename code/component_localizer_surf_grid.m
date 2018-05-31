function [comp_psc, condition_names, component_names] = ...
    component_localizer_surf_grid(  us, grid_roi, grid_spacing_mm, ...
    component_info, test_info, fwhm, varargin )

% 2016-09-09: Last modified, Sam NH
%
% 2017-01-09: Modified to accomodate combining data across runs in the first
% level analysis
% 
% 2017-01-27: Made pseudoinverse optional (i.e. to make it compatible with
% inverse analyses), added detrending default parameter, made 0 permutations the
% default
% 
% 2017-03-20: Very minor change: got rid of "remove_unspecified_trials"

% optional arguments
I.verbose = true;
I.anatomical_mask = '';
I.thresh_logP_residual_permtest = -inf;
I.pinv = true;
I.keyboard = false;
I = parse_optInputs_keyvalue(varargin, I);

if I.keyboard
    keyboard;
end

% default parameters
test_info = ...
    default_test_parameters(test_info, us);
component_info = ...
    default_component_parameters(component_info, us);

% add a permutation test in order to evaluate significance
if I.thresh_logP_residual_permtest > -inf
    component_info.n_perms = 100;
end

% re-do first level analyses if requested
if ~component_info.combine_runs_before_fla && component_info.overwrite_first_level
    for i = 1:length(component_info.runs)
        glm_surf_grid(...
            component_info.exp, us, component_info.runtype, ...
            fwhm, component_info.analysis_name, ...
            grid_spacing_mm, grid_roi, ...
            'n_perms', component_info.n_perms, ...
            'analysis_type', component_info.analysis_type, ...
            'runs', component_info.runs(i), 'plot_surf', false,...
            'plot_reliability', false, ...
            'overwrite_first_level', true, ...
            'whiten', component_info.whiten, ...
            'tsnr_threshold', component_info.tsnr_threshold);
    end
    component_info.overwrite_first_level = false;
end

for k = 1:length(test_info.runs) % loop through runs
    
    if I.verbose
        % print information about this localizer
        fprintf('Test: %s, run %d\n', test_info.runtype, test_info.runs(k));
        drawnow;
    end
    
    % matrix of psc values for each voxel
    % condition x voxel matrix
    [voxel_psc, condition_names] = psc_single_run(test_info, us, ...
        test_info.runs(k), fwhm, grid_spacing_mm, grid_roi);
    
    % create the mask
    n_voxels = size(voxel_psc,2);
    if isempty(I.anatomical_mask)
        mask = true(1,n_voxels);
    else
        mask = label2grid(I.anatomical_mask, grid_roi, grid_spacing_mm);
    end
    mask = mask > 0.99;
    
    % ensure non-independence
    if strcmp(test_info.exp, component_info.exp) ...
            && strcmp(test_info.runtype, component_info.runtype)
        localizer_runs_to_use = setdiff(...
            component_info.runs, test_info.runs(k));
    else
        localizer_runs_to_use = component_info.runs;
    end
        
    % component weights and significance values from a permutation test
    [comp_weights, logP_residual_permtest, component_names] = ...
        localizer_weights(component_info, us, localizer_runs_to_use, ...
        fwhm, grid_spacing_mm, grid_roi);
    
    % further select voxels based on significance values of permutation test
    if I.thresh_logP_residual_permtest > -inf;
        mask = mask & logP_residual_permtest > I.thresh_logP_residual_permtest;
    end
    
    if k == 1
        % component response matrix, initialization
        comp_psc = nan(length(test_info.runs), length(condition_names),...
            length(component_names));
    end
    
    % measure psc
    mask = mask & any(~isnan(voxel_psc),1) & all(~isnan(comp_weights),1);
    
    % apply weights
    if I.pinv 
        W = pinv(comp_weights(:,mask));
    else
        W = comp_weights(:,mask)';
    end
    
    try
        comp_psc(k,:,:) = voxel_psc(:,mask) * W;
        clear W;
    catch me
        print_error_message(me);
        keyboard;
    end
    
end

% helper function that find the appropriate file with p-values
function [beta_one_per_regressor, logP_residual_permtest, component_names] = ...
    localizer_weights(component_info, us, localizer_runs_to_use, ...
    fwhm, grid_spacing_mm, grid_roi) %#ok<STOUT>

[MAT_file_second_level, MAT_files_first_level, ...
    perm_MAT_file_second_level, perm_MAT_files_first_level, P] ...
    = glm_surf_grid(...
    component_info.exp, us, component_info.runtype, ...
    fwhm, component_info.analysis_name, ...
    grid_spacing_mm, grid_roi, ...
    'n_perms', component_info.n_perms, ...
    'analysis_type', component_info.analysis_type, ...
    'runs', localizer_runs_to_use, 'plot_surf', false,...
    'plot_reliability', false, ...
    'overwrite_first_level', component_info.overwrite_first_level, ...
    'overwrite_second_level', component_info.overwrite_second_level, ...
    'whiten', component_info.whiten, ...
    'combine_runs_before_fla', component_info.combine_runs_before_fla, ...
    'tsnr_threshold', component_info.tsnr_threshold);

use_first_level = (length(localizer_runs_to_use) == 1 ...
    || component_info.combine_runs_before_fla);

use_permtest = component_info.n_perms > 0;

if use_first_level
    fprintf('First level file\n%s\n', MAT_files_first_level{1}); drawnow;
    assert(length(MAT_files_first_level)==1);
    load(MAT_files_first_level{1}, 'beta_one_per_regressor');
    if use_permtest
        load(perm_MAT_files_first_level{1}, 'logP_residual_permtest');
    else
        logP_residual_permtest = [];
    end
else
    fprintf('Second level file\n%s\n', MAT_file_second_level); drawnow;
    load(MAT_file_second_level, 'beta_one_per_regressor');
    if use_permtest
        load(perm_MAT_file_second_level, 'logP_residual_permtest');
    else
        logP_residual_permtest = [];
    end
end

component_names = P.regressor_names;

function component_info = default_component_parameters(...
    component_info, us)

% runs to use
if ~isfield(component_info, 'runs')
    component_info.runs = read_runs(...
        component_info.exp, us, component_info.runtype);
end

% whether or not to overwrite first or second level files
if ~isfield(component_info, 'overwrite_first_level')
    component_info.overwrite_first_level = false;
end
if ~isfield(component_info, 'overwrite_second_level')
    component_info.overwrite_second_level = false;
end
if isfield(component_info, 'overwrite') && component_info.overwrite
    component_info.overwrite_first_level = true;
    component_info.overwrite_second_level = true;
end

% whether or not to combine across runs before performing first level analysis
if ~isfield(component_info, 'combine_runs_before_fla')
    component_info.combine_runs_before_fla = false;
end

% by default no permutations
if ~isfield(component_info, 'n_perms')
    component_info.n_perms = 0;
end

% by default, don't whiten
if ~isfield(component_info, 'whiten')
    component_info.whiten = false;
end

% by default, don't use a tsnr threshold
if ~isfield(component_info, 'tsnr_threshold')
    component_info.tsnr_threshold = 0;
end