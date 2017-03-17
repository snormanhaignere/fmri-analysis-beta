function G = template_surface_grid(subjid, grid_roi, grid_spacing_mm, varargin)

% Returns an empty grid file for a given subject. Grid voxels within the
% grid_roi have a value of 1. Grid voxels outside the grid_roi have a value of
% NaN.
% 
% 2017-03-10: Created, Sam NH

I.plot_figures = false;
I.overwrite = false;

global root_directory;
freesurfer_directory = [root_directory '/freesurfer/' subjid];

% flat patches
patch_rh = [freesurfer_directory '/surf/rh.cortex.patch.flat'];
patch_lh = [freesurfer_directory '/surf/lh.cortex.patch.flat'];

% roi files
roi_rh_label = [freesurfer_directory '/label/rh.' grid_roi '.label'];
roi_lh_label = [freesurfer_directory '/label/lh.' grid_roi '.label'];

% interpolates area map, but any map will do because we will discard the
% interpolated values later
surf_rh = [freesurfer_directory '/surf/rh.area.mgh'];
surf_lh = [freesurfer_directory '/surf/lh.area.mgh'];

% directory and file to save grid to
grid_directory = [freesurfer_directory '/grids'];
if ~exist(grid_directory, 'dir'); mkdir(grid_directory); end
grid_file = [...
    grid_directory '/grid-' num2str(grid_spacing_mm) 'mm_' grid_roi '.mat'];

if ~exist(grid_file, 'file') || I.overwrite
    
    % interpolate to the grid
    G = interp_from_surface_to_grid(...
        surf_rh, surf_lh, patch_rh, patch_lh, ...
        roi_rh_label, roi_lh_label, grid_spacing_mm, I.plot_figures);
    
    % erase area data 
    for q = 1:2
        G.grid_data{q}(~isnan(G.grid_data{q})) = 1;
    end
    
    % save
    save(grid_file, 'G');
    
else
    
    load(grid_file, 'G');
    
end


