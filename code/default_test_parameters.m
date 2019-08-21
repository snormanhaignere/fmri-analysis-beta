function test_info = default_test_parameters(test_info, us)

% Helper function fro roi_surf_grid.m and component_localizer_surf_grid.m
% 
% 2016-09-02: Separated out from roi_surf_grid.m, Sam NH
% 
% 2016-09-04: No longer defaults to older signal averaging analysis, condition
% names are not included in the defaults
% 
% 2019-08-10: Updated to include onset / offset delay

% runs to use
if ~isfield(test_info, 'runs')
    test_info.runs = read_runs(...
        test_info.exp, us, test_info.runtype);
end

% test info
if ~isfield(test_info, 'overwrite')
    test_info.overwrite = false;
end

% by default no tsnr threshold
if ~isfield(test_info, 'tsnr_threshold')
    test_info.tsnr_threshold = 0;
end

% default onsets and offsets for signal averaging
if ~isfield(test_info, 'onset_delay')
    test_info.onset_delay = 5;
end
if ~isfield(test_info, 'offset_delay')
    test_info.offset_delay = 1;
end

% number of time  points
if ~isfield(test_info, 'n_tps')
    test_info.n_tps = 1;
end