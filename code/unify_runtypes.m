function unified_grid_files = unify_runtypes(...
    exp, us, original_runtypes, unified_runtype, varargin)

% Combines surface grid files across sets of runs and for each group creates a
% new run of a new combined runtype.
%
% 2017-02-23: Created, Sam NH

global root_directory;

% optional input arguments
I.overwrite = false;
I.keyboard = false;
I.original_para_prefix = original_runtypes;
I.unified_para_prefix = unified_runtype;
I.surf_grid_params = {};

% by default include all runs
I.original_runs = cell(size(original_runtypes));
for i = 1:length(original_runtypes)
    I.original_runs{i} = read_runs(exp, us, original_runtypes{i});
end

% by default runs of the unified run type are given sequential numbers
% i.e. {[1,2,3], [4,5], [7,8,9], ...}
I.unified_runs = cell(size(I.original_runs));
n_total_runs = 0;
for i = 1:length(I.original_runs)
    I.unified_runs{i} = n_total_runs + (1:length(I.original_runs{i}));
    n_total_runs = n_total_runs + length(I.original_runs{i});
end

% reset optional inputs
I = parse_optInputs_keyvalue(varargin, I);

if I.keyboard
    keyboard;
end

%% Combine surface files

% run order file with the new orders
run_order_file = ...
    [root_directory '/' exp '/data/brain/runorders/' ...
    unified_runtype '_us' num2str(us) '.txt'];
if ~exist(run_order_file, 'file') || I.overwrite
    fid = fopen(run_order_file, 'w');
else
    fid = NaN;
end

unified_grid_files = cell(1,length(I.original_runs));
for i = 1:length(I.original_runs)
    
    [runs, seqid, scanid] = read_runs(exp, us, original_runtypes{i});
    for j = 1:length(I.original_runs{i})
        
        % add a line to run_order_file
        if ~isnan(fid)
            xi = I.original_runs{i}(j) == runs;
            assert(sum(xi)==1);
            fprintf(fid, '%d %d %d\n', ...
                I.unified_runs{i}(j), seqid(xi), scanid(xi));
        end
        
        % soft link to para file
        original_para_file = ...
            [root_directory '/' exp '/data/para/usub' num2str(us) ...
            '/' I.original_para_prefix{i}  ...
            '_r' num2str(I.original_runs{i}(j)) '.par'];
        unified_para_file = ...
            [root_directory '/' exp '/data/para/usub' num2str(us) ...
            '/' I.unified_para_prefix  ...
            '_r' num2str(I.unified_runs{i}(j)) '.par'];
        if ~exist(unified_para_file, 'file') || I.overwrite
            unix(['ln -s -f ' original_para_file ' ' unified_para_file]);
        end
        
        % soft link to scan parameter file
        % slightly more complicated because of need to determine if there is a
        % run-specific parameter file
        scan_parameter_directory = ...
            [root_directory '/' exp '/data/brain/scan-parameters'];
        original_scan_parameter_file = ...
            [scan_parameter_directory '/' original_runtypes{i} ...
            '_us' num2str(us) '_r' num2str(I.original_runs{i}(j)) '.txt'];
        if ~exist(original_scan_parameter_file, 'file')
            original_scan_parameter_file = ...
                [scan_parameter_directory '/' original_runtypes{i} '.txt'];
        end
        unified_scan_parameter_file = ...
            [scan_parameter_directory '/' unified_runtype ...
            '_us' num2str(us) '_r' num2str(I.unified_runs{i}(j)) '.txt'];
        if ~exist(unified_scan_parameter_file, 'file') || I.overwrite
            unix(['ln -s -f ' original_scan_parameter_file ...
                ' ' unified_scan_parameter_file]);
        end
        
        % link preprocessing directory
        original_preprocessing_directory = ...
            [root_directory '/' exp '/analysis/preprocess/usub' num2str(us) ...
            '/' original_runtypes{i} '_r' num2str(I.original_runs{i}(j))];
        if ~exist(original_preprocessing_directory, 'dir')
            mkdir(original_preprocessing_directory)
        end
        unified_preprocessing_directory = ...
            [root_directory '/' exp '/analysis/preprocess/usub' num2str(us) ...
            '/' unified_runtype '_r' num2str(I.unified_runs{i}(j))];
        if ~exist(unified_preprocessing_directory, 'file') || I.overwrite
            unix(['ln -f -s ' original_preprocessing_directory ...
                ' ' unified_preprocessing_directory]);
        end
        
        % nifti data files
        original_nifti_file = ...
            [root_directory '/' exp '/data/brain/nifti/usub' num2str(us) ...
            '/' original_runtypes{i} '_r' num2str(I.original_runs{i}(j)) '.nii.gz'];
        unified_nifti_file = ...
            [root_directory '/' exp '/data/brain/nifti/usub' num2str(us) ...
            '/' unified_runtype '_r' num2str(I.unified_runs{i}(j)) '.nii.gz'];
        if exist(original_nifti_file, 'file') && ...
                (~exist(unified_para_file, 'file') || I.overwrite)
            unix(['ln -s -f ' original_nifti_file ' ' unified_nifti_file]);
        end
        
    end
end

