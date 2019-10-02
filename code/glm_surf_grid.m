function [MAT_file_second_level, MAT_files_first_level, ...
    perm_MAT_file_second_level, perm_MAT_files_first_level, ...
    P, analysis_directory, figure_directory, template_grid_file] = ...
    glm_surf_grid(exp, us, runtype, fwhm, analysis_name, ...
    grid_spacing_mm, grid_roi, varargin)

% 2016-08-27: Modified how optional arguments are handled, Sam NH
%
% 2016-08-31: Made the prefix of the para files an optional argument, Sam NH
%
% 2016-09-09: n_perms made an optional argument, permutation tests are saved as
% a separate MAT file, requiring small changes to this wrapper function, Sam NH
%
% 2016-09-09: This function now directly calls first level analysis script
% glm_event_regression.m, and handles whether or not to overwrite files, Sam NH
%
% 2016-09-09: Generalized glm_surf_grid so that it can also handled GLMs run on
% signal averaged responses, Sam NH
%
% 2016-11-18: Returns template_grid_file, Sam NH
%
% 2017-01-09: Made it possible to combine data across runs before performing a
% first level analysis, Sam NH
%
% 2017-01-27: Fixed up some issues with run combination. Made it possible to run
% inverse glm analyses, regressing voxels against a predicted response.
%
% 2017-02-24: Flag 'remove_unspecified_conditions' removed
% 
% 2019-08-20: Added ability to measure multiple timepoints in between onset
% and offset so as to measure a signal averaged timecourse (only for
% sigav-glm)

global root_directory;

%% Optional arguments

% optional arguments and defaults
I.overwrite = false;
I.overwrite_second_level = false;
I.overwrite_first_level = false;
I.plot_surf = false;
I.stat_to_plot = 'logP_permtest';
I.color_range = NaN;
I.plot_reliability = false;
I.runs = read_runs(exp, us, runtype);
I.n_perms = 0;
I.analysis_type = 'glm';
I.onset_delay = 5; % only applicable for signal averaging
I.offset_delay = 1; % only applicable for signal averaging
I.whiten = false;
I.renderer = 'opengl';
I.combine_runs_before_fla = false;
I.keyboard = false;
I.tsnr_threshold = 0;
I.grid_flag = '';
I.write_surf_file = false;
I.n_tps = 1;
[I, C] = parse_optInputs_keyvalue(varargin, I, 'empty_means_unspecified', true);
if I.overwrite
    I.overwrite_first_level = true;
    I.overwrite_second_level = true;
end

if I.keyboard
    keyboard;
end

%% String identifying the parameters of this analysis

param_idstring = ['fsaverage_smooth-' num2str(fwhm) 'mm' ...
    '_' 'grid-' num2str(grid_spacing_mm) 'mm_' grid_roi];

if ~isempty(strfind(I.analysis_type, 'sigav')) && C.onset_delay
    param_idstring = [param_idstring '_ons' num2str(I.onset_delay)];
end

if ~isempty(strfind(I.analysis_type, 'sigav')) && C.offset_delay
    param_idstring = [param_idstring '_off' num2str(I.offset_delay)];
end

if C.tsnr_threshold
    param_idstring = [param_idstring '_tsnr' num2str(I.tsnr_threshold)];
end

if ~isempty(I.grid_flag)
    param_idstring = [param_idstring '_' I.grid_flag];
end

if C.n_tps
    param_idstring = [param_idstring '_ntps' num2str(I.n_tps)];
end

%% Directories / setup

% analysis directory
analysis_directory = [root_directory  '/' exp '/analysis' ...
    '/' I.analysis_type '/' analysis_name ...
    '/' param_idstring '/usub' num2str(us)];

% condition weight file
parameter_file = [root_directory '/' exp '/analysis/' I.analysis_type ...
    '/' analysis_name '/parameters.mat'];
P = load(parameter_file);
P = glm_default_parameters(P);

% analysis directory
figure_directory = strrep(analysis_directory, 'analysis', 'figures');

% create this directories if not present
if ~exist(analysis_directory, 'dir')
    mkdir(analysis_directory);
end
if ~exist(figure_directory, 'dir')
    mkdir(figure_directory);
end

% prefix for the paradigm file
if isfield(P, 'para_prefix')
    para_prefix = P.para_prefix;
else
    para_prefix = runtype;
end

%% Ensure all surface files exist

for i = 1:length(I.runs)
    downsample_surface_timecourses(exp, us, runtype, I.runs(i), fwhm, ...
        grid_spacing_mm, grid_roi, 'plot', false);
end

%% Combine across runs first

if I.combine_runs_before_fla
    
    % combining data across runs
    combined_runtype = [runtype '_combined_detrend' num2str(P.detrend)];
    combined_para_prefix = [para_prefix '_combined_detrend' num2str(P.detrend)];
    combine_runs_wrapper(...
        exp, us, runtype, {I.runs}, fwhm, grid_roi, grid_spacing_mm, ...
        'detrend', P.detrend, 'original_para_prefix', para_prefix, ...
        'combined_runtype', combined_runtype, ...
        'combined_para_prefix', combined_para_prefix, ...
        'overwrite', I.overwrite_first_level, ...
        'list_all_runs_in_fname', true, 'create_run_order_file', false);
    
    % no need to detrend further
    % creates a new parameter file with detrend=0
    % file deleted at the end of this script
    P.detrend = 0;
    [a,b,c] = fileparts(parameter_file);
    parameter_file = [a '/' b '_0detrend' c];
    save(parameter_file, '-struct', 'P');
    delete_parameter_file = true;
    
    % sets of runs for the first level analysis
    fla_run_sets = {I.runs};
    
    % set run type and para prefix to their new versions
    runtype = combined_runtype;
    para_prefix = combined_para_prefix;
    
else
    
    fla_run_sets = cell(1, length(I.runs));
    for i = 1:length(I.runs)
        fla_run_sets{i} = I.runs(i);
    end
    delete_parameter_file = false;
    
end

%% First level analysis separately for each run

% create cell struction with para files and data files
n_run_sets = length(fla_run_sets);
para_files = cell(1, n_run_sets);
data_matrix_files = cell(1, n_run_sets);
MAT_files_first_level = cell(1, n_run_sets);
perm_MAT_files_first_level = cell(1, n_run_sets);

for i = 1:n_run_sets
    
    % string with all of the runs to be analyzed in this set
    run_idstring = sprintf('%d',fla_run_sets{i});
    
    % print string to command window
    fprintf('First level analysis of run(s) %s\n', run_idstring); drawnow;
    
    % para file
    para_files{i} = [root_directory '/' exp '/data/para/usub' num2str(us) ...
        '/' para_prefix  '_r' run_idstring '.par'];
    
    % TR
    TR = read_functional_scan_parameters(exp,us,runtype,fla_run_sets{i}); %#ok<NASGU>
    
    % preprocessing directory with files in fsaverage space
    preproc_fsaverage_directory = [root_directory '/' exp '/analysis' ...
        '/preprocess/usub' num2str(us) '/' ...
        runtype '_r' run_idstring '/myfsaverage'];
    
    % input surface grid
    grid_file = [preproc_fsaverage_directory '/' ...
        'smooth-' num2str(fwhm) 'mm' '_' ...
        'grid-' num2str(grid_spacing_mm) 'mm' '_' grid_roi '.mat'];
    if ~isempty(I.grid_flag)
        grid_file = strrep(grid_file, '.mat', ['_' I.grid_flag '.mat']);
    end
    
    % reformated data matrix to use as input to the GLM analysis below
    data_matrix_files{i} = ...
        strrep(grid_file, '.mat', '_unwrapped_data_matrix.mat');
    
    % a blank template of the surface grid for future reference
    template_grid_file = [analysis_directory '/template_grid.mat'];
    
    % create the data matrix and blank template, if not present already
    if ~exist(template_grid_file, 'file') ...
            || ~exist(data_matrix_files{i}, 'file') || I.overwrite_first_level
        
        fprintf('Converting grid files to data matrix...\n'); drawnow;
        
        % reformat to voxel x datapoint/timepoint
        load(grid_file, 'G');
        data_matrix = grid2matrix(G); %#ok<NASGU>
        
        % save
        save(data_matrix_files{i}, 'data_matrix', 'TR', 'G', '-v7.3');
        
        % set all values to NaN and remove third dimension
        G.grid_data{1} = nan(size(G.grid_data{1},1), size(G.grid_data{1},2));
        G.grid_data{2} = nan(size(G.grid_data{2},1), size(G.grid_data{2},2));
        save(template_grid_file, 'G', '-v7.3');
        
    end
    
    % file to save results of individual run analysis
    MAT_files_first_level{i} = ...
        [analysis_directory '/' runtype '_r' run_idstring '.mat'];
    
    % first level MAT files with permuted stats
    if I.n_perms > 0
        perm_MAT_files_first_level{i} = ...
            [analysis_directory '/' runtype '_r' run_idstring ...
            '_' num2str(I.n_perms) 'perms.mat'];
    end
    
    % if multiple runs, indicate this is a first level analysis
    if length(fla_run_sets{i})>1
        MAT_files_first_level{i} = ...
            strrep(MAT_files_first_level{i}, '.mat', '_fla.mat');
        if I.n_perms > 0
            perm_MAT_files_first_level{i} = ...
                strrep(perm_MAT_files_first_level{i}, '.mat', '_fla.mat');
        end
    end
    
    % check if output files already exist
    if ~exist(MAT_files_first_level{i}, 'file') ...
            || (I.n_perms > 0 && ~exist(perm_MAT_files_first_level{i}, 'file'))...
            || I.overwrite_first_level
        
        switch I.analysis_type
            case 'glm'
                
                fprintf('GLM...\n'); drawnow;
                assert(I.n_tps == 1);
                
                % add white matter regressors
                X_nuissance = [];
                if P.n_whitematter_PCs > 0
                    error('Need to set this up to work with run sets');
                    PCs = whitematter_PCs(exp, us, runtype, r, 'motcorr', 'bbreg');
                    X_nuissance = [X_nuissance; PCs(:,1:P.n_whitematter_PCs)];
                end
                
                % save nuissance regressors
                if ~isempty(X_nuissance)
                    nuissance_regressor_file = [analysis_directory '/' runtype ...
                        '_r' run_idstring '_nuissance.mat'];
                    save(nuissance_regressor_file, 'X_nuissance');
                else
                    nuissance_regressor_file = [];
                end
                
                % first level regression
                glm_event_regression(data_matrix_files{i}, para_files{i}, ...
                    parameter_file, MAT_files_first_level{i}, ...
                    'n_perms', I.n_perms, 'nuissance_regressor_file', ...
                    nuissance_regressor_file, 'tsnr_threshold', I.tsnr_threshold);
                
            case 'sigav-glm'
                
                fprintf('Signal averaging/GLM...\n'); drawnow;
                
                % first level regression
                sigav_glm(data_matrix_files{i}, para_files{i}, ...
                    parameter_file, MAT_files_first_level{i}, ...
                    'onset_delay', I.onset_delay, 'offset_delay', I.offset_delay,...
                    'n_perms', I.n_perms, 'whiten', I.whiten, ...
                    'tsnr_threshold', I.tsnr_threshold, ...
                    'n_tps', I.n_tps);
                
            case 'sigav-inverse-glm'
                
                fprintf('Signal averaging/Inverse GLM...\n'); drawnow;
                assert(I.n_tps == 1);

                % first level regression
                sigav_inverse_glm(data_matrix_files{i}, para_files{i}, ...
                    parameter_file, MAT_files_first_level{i}, ...
                    'onset_delay', I.onset_delay, 'offset_delay', I.offset_delay,...
                    'n_perms', I.n_perms, 'tsnr_threshold', I.tsnr_threshold);
                
            otherwise
                error('No matching case for analysis type "%s"\n', I.analysis_type);
                
        end
    end
end

%% Second level analysis

if length(I.runs) > 1 && ~I.combine_runs_before_fla ...
        && ~strcmp(I.analysis_type, 'sigav-inverse-glm')
    
    % file to save results of second level analysis pooling across runs
    MAT_file_second_level = [analysis_directory '/' runtype '_r' sprintf('%d', I.runs) '.mat'];
    
    % perform second level analysis, first check if output already exists
    if ~exist(MAT_file_second_level, 'file') || I.overwrite_second_level
        second_level_ols(MAT_files_first_level, MAT_file_second_level);
    end
    
    % second level MAT files with permuted stats
    if I.n_perms > 0
        perm_MAT_file_second_level = [analysis_directory ...
            '/' runtype '_r' sprintf('%d', I.runs) '_' num2str(I.n_perms) 'perms.mat'];
        
        % perform the second level analysis
        if ~exist(perm_MAT_file_second_level, 'file') || I.overwrite_second_level
            second_level_permtest(MAT_files_first_level, ...
                perm_MAT_files_first_level, perm_MAT_file_second_level);
        end
    else
        perm_MAT_file_second_level = [];
    end
    
else
    
    MAT_file_second_level = [];
    perm_MAT_file_second_level = [];
    
end

%% Analysis concatenating across runs

% for i = 1:length(I.runs)
%
%                 % first level regression
%                 sigav_glm(data_matrix_files{i}, para_files{i}, ...
%                     parameter_file, MAT_files_first_level{i}, ...
%                     'onset_delay', I.onset_delay, 'offset_delay', I.offset_delay,...
%                     'n_perms', I.n_perms, 'whiten', I.whiten, ...
%                     'remove_unspecified_trials', I.remove_unspecified_trials);

%% Reliability measures

if length(MAT_files_first_level) > 1 && I.plot_reliability
    
    % plot reliability of contrast across runs
    glm_contrast_map_reliability(MAT_files_first_level,...
        analysis_directory, figure_directory, runtype, 'overwrite', I.overwrite);
    
    % reliability of individual voxel responses across regressors
    glm_voxel_reliability(MAT_files_first_level,...
        analysis_directory, figure_directory, runtype, 'overwrite', I.overwrite);
    
    % % plot reliability of response profile across runs
    % glm_regressor_response_reliability(MAT_files_first_level,...
    %     analysis_directory, figure_directory, 'overwrite', I.overwrite);
    
end

%% Plot surface

if delete_parameter_file
    delete(parameter_file);
end

if ~I.plot_surf
    return;
end

assert(I.n_tps == 1);

% select first or second level analysis
if ~isempty(MAT_file_second_level)
    matfile = MAT_file_second_level;
else
    matfile = MAT_files_first_level{1};
end

% stats from permutation test are stored in a separate file
if strfind(I.stat_to_plot, 'permtest');
    [parent_directory, MAT_file_noext] = fileparts(matfile);
    matfile = [parent_directory '/' ...
        MAT_file_noext '_' num2str(I.n_perms) 'perms.mat'];
    clear parent_directory matfile_noext;
end

% load the stat
X = load(matfile, I.stat_to_plot);

% shape to grid
load(template_grid_file, 'G');
G = matrix2grid(X.(I.stat_to_plot)', G);
surf = grid2surface(G);

% color range to plot
if ~isnan(I.color_range)
    color_range = I.color_range;
    assert(length(color_range)==2);
else
    switch I.stat_to_plot
        case {'logP_ols', 'logP_permtest'};
            color_range = [-5 5];
        case {'beta_contrast', 'beta_one_per_regressor'};
            color_range = [-3 3];
        case {'logP_residual_permtest'}
            color_range = [0 5];
        otherwise
            error('No matching case');
    end
end

n_maps = size(surf,3);
for i = 1:n_maps
    hemis = {'rh','lh'};
    for q = 1:2
        
        if ~exist('figh', 'var') % Plot figures with outlines of standard ROIs
            close all;
            figh = figure;
            pos = get(figh,'Position');
            set(figh, 'Renderer', I.renderer);
            set(figh, 'Position', [pos(1:2), 800 800]);
        end
        
        % plot surface map
        plot_fsaverage_1D_overlay(surf(:,q,i),hemis{q},'parula',color_range,figh);
        
        try
            % save to file
            switch I.stat_to_plot
                case {'logP_ols', 'beta_contrast'}
                    fname_substring = [P.contrast_names{i}];
                case {'logP_permtest'}
                    fname_substring = [P.contrast_names{i} '_' num2str(I.n_perms) 'perms'];
                case {'beta_one_per_regressor'}
                    fname_substring = [P.regressor_names{i}];
                case {'logP_residual_permtest'}
                    fname_substring = [num2str(I.n_perms) 'perms'];
                otherwise
                    error('No matching case');
            end
            figure_file = [figure_directory '/' runtype '_r' sprintf('%d', I.runs) '_pmap_' I.stat_to_plot '_' ...
                fname_substring '_' hemis{q} '_colrange_' ...
                num2str(color_range(1)) '_' num2str(color_range(2)) '.png'];
            export_fig(figure_file,'-png','-r100','-nocrop');
            
            if I.write_surf_file
                surf_file = [figure_directory '/' runtype '_r' sprintf('%d', I.runs) '_pmap_' I.stat_to_plot '_' ...
                    fname_substring '_' hemis{q} '.mgh'];
                MRIwrite_surface(surf(:,q,i), surf_file, hemis{q});
            end
            
        catch me
            print_error_message(me);
            keyboard;
        end
    end
end