function [color_data, patch_handle, light_handle] = ...
    plot_fsaverage_1D_overlay_v2(surface_values, hemi, varargin)

% function [color_data, patch_handle, light_handle] = plot_fsaverage_1D_overlay_v2(surface_values, hemi, color_map_to_plot, color_range, figh, varargin)
% 
% Plots a surface overlay on the fsaverage template brain using the matlab
% function patch. Version 2 has optional arguments specified as key-value
% pairs.
% 
% -- Inputs --  
% 
% surface_values: vector of values, one per vertex, to plot, NaN values are not plotted
% 
% hemi: whether the left or right hemisphere is being plotted
% 
% color_map_to_plot (optional): the name of the colormap to use (default is 'parula' if not specified), can
% also be a N x 3 matrix of values to interpolate within the color range.
% 
% color_range (optional): range of values to plot, [lower_bound, upperbound] (if not specified the central 95 of the distribution of values is plotted)
% 
% figh (optional): matlab handle of the figure to plot the surface in, if unspecified a new figure
% handle is created (i.e. figh = figure)
% 
% -- Outputs --
% 
% color_data: N x 3 matrix specifying the RGB color of each vertex
% 
% patch_handle: handle to the patch object created
% 
% light_handle: handle to the light object created
% 
% -- Example: Plots significance map for music component discovered by ICA --  
% 
% hemi = 'rh';
% surf = MRIread(['/mindhive/nklab/u/svnh/fmri-analysis/test_data/' hemi '.ICA_pmap_music.mgz']);
% surface_values = surf.vol;
% colormapname = 'parula';
% plot_fsaverage_1D_overlay_v2(surface_values, hemi, 'colormap', colormapname);
% 
% Modified by Sam NH on 9/3/2015

global root_directory;

clear I;
I.colormap = 'parula';
I.color_range = [];
I.logscale = false;
I.plot = true;
I.figh = matlab.ui.Figure.empty;
I = parse_optInputs_keyvalue(varargin, I, 'ignore_mismatch_class', {'colormap'});

if isempty(I.color_range)
    [Nx,x] = hist(surface_values,100);
    Cx = cumsum(Nx/sum(Nx));
    [~,xi] = unique(Cx);
    x = x(xi);
    Cx = Cx(xi);
    color_range = interp1(Cx,x,[0.025 0.975]);
else
    color_range = I.color_range;
end

if I.logscale
    surface_values = log2(surface_values);
    color_range = log2(color_range);
end

% read vertices and faces
nvertices = 163842; 
[vertices, faces] = freesurfer_read_surf([root_directory '/' 'freesurfer/myfsaverage/surf/' hemi '.inflated'], false);

% set default color data based on gyral/sulcal divisions
curv = read_curv([root_directory '/' 'freesurfer/myfsaverage/surf/' hemi '.curv']);
color_data = 0.5*ones(nvertices,3);
color_data(curv>0,:) = 0.3;

% read in a colormap
if ischar(I.colormap)
    h = figure;
    cmap = colormap(I.colormap);
    close(h);
else
    cmap = I.colormap;
end

% interpolate surface values to colormap
surface_values_bounded = surface_values;
surface_values_bounded(surface_values_bounded < color_range(1)) = color_range(1);
surface_values_bounded(surface_values_bounded > color_range(2)) = color_range(2);
x = linspace(color_range(1),color_range(2),size(cmap,1))';
for i = 1:3
    color_data(~isnan(surface_values_bounded),i) = interp1(x, cmap(:,i), surface_values_bounded(~isnan(surface_values_bounded))','pchip');
end

% return after calculating color data without plotting
if ~I.plot
    return;
end

% create figure of specified size
if nargin < 5 || isempty(I.figh)
    I.figh = figure;
    pos = get(I.figh,'Position');
    set(I.figh, 'Position', [pos(1:2), 800 800]);
    clf(I.figh);
else
    clf(I.figh);
end

% create the patch object
patch_handle = patch('vertices', vertices, 'Faces', faces, 'FaceVertexCData', color_data,'FaceLighting','gouraud','SpecularStrength',0,'DiffuseStrength',0.7);
shading interp;

% adjust viewing angle 
switch hemi
    case {'rh'}
%         view([115 10])
        camup([0 0.5 1]);
        camva(4.532)
        campos([1.3604e+03 1.0309e+03 363.9286]);
        camtarget([0 15 -10]);
        xlim(2.2*40*[-1 1]); ylim(2.2*65*[-1 1]); zlim(2.2*55*[-1 1]);
        light_handle = camlight('right','infinite');
        set(light_handle, 'Position', [0.33 0.33 0.33]);

    case {'lh'}
%         view([-105 0]);
        camup([0 0.5 1]);
        % camzoom(2.2);
        camtarget([-10 10 -5]);
        %         campos([-1449.9 631.312 363.929]);
        campos([-1449.9*1 631.312*1.3 363.929*1]);
        camva(4);
        % campos([1.3604e3 1.0309e3 0.3639e3])
        xlim(2.2*40*[-1 1]); ylim(2.2*65*[-1 1]); zlim(2.2*55*[-1 1]);
        light_handle = camlight('left','infinite');
        set(light_handle, 'Position', [-1, 1, 0.33]);
        
    otherwise
        error('hemi should be "rh" or "lh", not %s',hemi);
end

% colorbar
colormap(cmap);
colorbar_handle = colorbar('Location','South');
colorbar_labels_str = cell(1,5);
colorbar_labels_num = linspace(color_range(1),color_range(2),5);
if I.logscale
    colorbar_labels_num = 2.^colorbar_labels_num;
end
for i = 1:5
    num_sig_digits = round(log10(colorbar_labels_num(i)));
    if num_sig_digits > 3 || num_sig_digits < -3
        colorbar_labels_str{i} = sprintf('%.2e',colorbar_labels_num(i));
    else
        colorbar_labels_str{i} = sprintf('%.2f',colorbar_labels_num(i));
    end
end
ca = caxis;
set(colorbar_handle, 'XTick', linspace(ca(1),ca(2),5), 'XTickLabel', colorbar_labels_str,'FontSize',20,'Position',[0.1469 0.05 0.7438 0.0312]);
