function global_mean_center(exp, us, runtype, r, fwhm, ...
    grid_spacing_mm, grid_roi, fwhm_global_mean, varargin)

global root_directory;
project_directory = [root_directory '/' exp];

I.overwrite = false;
I.keyboard = false;
I.plot_smoothed_map = false;
I = parse_optInputs_keyvalue(varargin, I);

% input/output files
analysis_directory = [project_directory '/analysis/preprocess'...
    '/usub' num2str(us) '/' runtype '_r' num2str(r) '/myfsaverage'];
fname_noext = [analysis_directory ...
    '/smooth-' num2str(fwhm) 'mm_grid-' ...
    num2str(grid_spacing_mm) 'mm_' grid_roi];
input_file = [fname_noext '.mat'];
smooth_file = [fname_noext '_gm-' num2str(fwhm_global_mean) 'mm.mat'];
output_file = [fname_noext '_gmc-' num2str(fwhm_global_mean) 'mm.mat'];

if ~exist(smooth_file, 'file') || ~exist(output_file, 'file') || I.overwrite
    
    % load data
    load(input_file, 'G');
    G_orig = G;
    clear G;
    
    % enter debug mode
    if I.keyboard
        keyboard;
    end
    
    % compute the global mean
    G_smooth = G_orig;
    for i = 1:2
        G_smooth.grid_data{i} = mysmooth2(...
            G_orig.grid_data{i}, grid_spacing_mm, fwhm_global_mean, ...
            'plot_effect', double(I.plot_smoothed_map));
    end
    G = G_smooth;
    save(smooth_file, 'G');
    clear G;

    % center using global mean
    G_center = G_orig;
    for i = 1:2
        mean_of_timecourse = mean(G_smooth.grid_data{i},3);
        G_center.grid_data{i} = bsxfun(@plus, G_orig.grid_data{i} - G_smooth.grid_data{i}, mean_of_timecourse);
        clear mean_of_timecourse;
    end
    G = G_center;
    save(output_file, 'G');
    clear G;
    
end