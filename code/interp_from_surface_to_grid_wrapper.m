function G = interp_from_surface_to_grid_wrapper(surf_rh, surf_lh, grid_roi, ...
    grid_spacing_mm, subjid, plot_figures)

% Wrapper for interp_from_surface_to_grid.m that makes it easier to resample a
% surface to a grid.
% 
% 2016-01-24: Created, Sam NH

if nargin < 4 || isempty(grid_spacing_mm)
    grid_spacing_mm = 1.5;
end

if nargin < 5 || isempty(subjid)
    subjid = 'myfsaverage';
end

if nargin < 6 || isempty(plot_figures)
    plot_figures = false;
end

global root_directory;
freesurfer_directory = [root_directory '/freesurfer/' subjid];

% flat patches
patch_rh = [freesurfer_directory '/surf/rh.cortex.patch.flat'];
patch_lh = [freesurfer_directory '/surf/lh.cortex.patch.flat'];

% roi files
roi_rh_label = [freesurfer_directory '/label/rh.' grid_roi '.label'];
roi_lh_label = [freesurfer_directory '/label/lh.' grid_roi '.label'];

% interpolate to the grid
G = interp_from_surface_to_grid(...
    surf_rh, surf_lh, patch_rh, patch_lh, ...
    roi_rh_label, roi_lh_label, grid_spacing_mm, plot_figures);
