function [mat_file_with_all_useful_statistics, p_rh_inflated_cluster_thresh_file, p_lh_inflated_cluster_thresh_file] = tla_permutation_clustcorr(P,volume_or_surface,tla_directory_name,varargin)

% function mat_file_with_all_useful_statistics = tla_permutation_clustcorr(P,volume_or_surface,tla_directory_name,varargin)
% 
% Cluster-corrects voxelwise statistics calculated by tla_permutation_voxelwise.
% A distribution of cluster sizes is calculated using the permuted orders, analagous
% to standard cluster-correction based on the family-wise error.
% 
% -- Example from Amusia experiment --
% clear P;
% usubs_amusia  = [45,49,51,53,55,57,59,71,73,75,171]';
% tla_directory_name = 'amusia_group_harm_vs_noise_5mm';
% for i = 1:11;
%     P(i).exp = 'amusia'; %#ok<*SAGROW>
%     P(i).us = usubs_amusia(i);
%     P(i).runtype = 'localizer';
%     P(i).runs = 1;
%     P(i).contrast = 'harm_vs_noise';
%     % P(i).lower_level_directory_name = 'smooth300mm_grid_hand-stp-stg_1.5mm';
%     P(i).lower_level_directory_name = 'smooth500mm_grid_hand-stp-stg_1.5mm_10whitematterPCs';
% end
% volume_or_surface = 'downsampled_surface';
% tla_permutation_voxelwise_stats(P,volume_or_surface,tla_directory_name,'n_smps',10e3);
% tla_permutation_clustcorr(P,volume_or_surface,tla_directory_name,'n_smps',10e3,'sigmap');
% 
% Last modified by Sam Norman-Haignere on 15-06-18

fprintf('v6'); drawnow;

% scripts directories
source_directory = strrep(which('fla_matlab.m'),'fla_matlab.m','');
addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
addpath([source_directory 'export_fig']);

% default number of samples used to compute p-values
n_smps = 10e3;
if optInputs(varargin, 'n_smps')
    n_smps = varargin{optInputs(varargin, 'n_smps')+1};
end

% p-value threshold, which determines the cluster size
voxel_pthresh = 3;
if optInputs(varargin, 'voxel_pthresh')
    voxel_pthresh = varargin{optInputs(varargin, 'voxel_pthresh') + 1};
end

% cluster-corrected p-value
cluster_pthresh = 1.3;
if optInputs(varargin, 'cluster_pthresh')
    cluster_pthresh = varargin{optInputs(varargin, 'cluster_pthresh') + 1};
end


% files
switch volume_or_surface
    case 'volume'
        
        error('Need to setup volume analysis');
        
    case 'surface'
        
        error('Need to setup surface analysis');
        
    case 'downsampled_surface'
        
        % freesurfer subject id
        if optInputs(varargin, 'monkey')
            subjid = [P(1).exp '_us' num2str(P(1).us)];
            tla_directory = [params('rootdir') 'freesurfer/' subjid      '/tla_matlab/' tla_directory_name '_downsampled_hash' DataHash(P) '/'];
        else
            tla_directory = [params('rootdir') 'freesurfer/' 'fsaverage' '/tla_matlab/' tla_directory_name '_downsampled_hash' DataHash(P) '/'];
        end
        
        if ~exist(tla_directory,'dir');
            mkdir(tla_directory);
        end

        % p-values based on counting the number of times the null
        % distribution exceeds the measured value
        voxelwise_p_file =  [tla_directory  'pstat_permutation_' num2str(n_smps) '.mat'];
        
        % p-values based on gaussian fits to the null distribution
        voxelwise_p_gaussfit_file =  [tla_directory  'pstat_permutation_' num2str(n_smps) '_gaussfit.mat'];
        voxelwise_p_gaussfit_rh_inflated_file =  [tla_directory  'rh.pstat_permutation_' num2str(n_smps) '_gaussfit.mgz'];
        voxelwise_p_gaussfit_lh_inflated_file =  [tla_directory  'lh.pstat_permutation_' num2str(n_smps) '_gaussfit.mgz'];
        
        mat_file_with_all_useful_statistics =  [tla_directory  'all_statistics_' num2str(n_smps) '.mat'];
        
        p_cluster_thresh_file =  [tla_directory  'pstat_permutation_' num2str(n_smps) '_cluster_corrected_voxP' num2str(voxel_pthresh) '_clustP' num2str(cluster_pthresh)  '.mat'];
        p_rh_inflated_cluster_thresh_file =  [tla_directory  'rh.pstat_permutation_' num2str(n_smps) '_cluster_corrected_voxP' num2str(voxel_pthresh) '_clustP' num2str(cluster_pthresh)  '.mgz'];
        p_lh_inflated_cluster_thresh_file =  [tla_directory  'lh.pstat_permutation_' num2str(n_smps) '_cluster_corrected_voxP' num2str(voxel_pthresh) '_clustP' num2str(cluster_pthresh)  '.mgz'];
        
    otherwise
        error('Error in tla_permutation: volume_or_surface flag must be either "volume", "surface", or "downsampled_surface"...');
        
end

if  ~exist(p_cluster_thresh_file, 'file') || optInputs(varargin, 'overwrite')
        
    load(mat_file_with_all_useful_statistics, 'zstat', 'pstat', 'pstat_gaussfit', 'zstat_gaussfit', 'nullmean','nullstd', 'permuted_lower_level_z_with_noNaNs', 'voxels_with_no_NaNs', 'nullmean_voxels_with_no_NaNs', 'nullstd_voxels_with_no_NaNs', 'smps');
    load(voxelwise_p_file);

    % initialize variables
    n_subjects = size(permuted_lower_level_z_with_noNaNs,2); %#ok<NODEF>
    n_voxels = length(zstat);
    sampled_lower_level_zmap = nan(n_subjects, sum(voxels_with_no_NaNs));
    sampled_tla_pmap_gaussfit = nan(1, n_voxels);
    rh_grid = nan(size(G.grid_data{1}));
    lh_grid = nan(size(G.grid_data{2}));
    n_clusters = nan(1,n_smps);
    max_cluster_size = zeros(1,n_smps);
    
    fprintf('Cluster statistics for permuted z-maps\n');
    tic;
    for i = 1:n_smps
        
        % z-map for this set of samples
        for j = 1:n_subjects
            % select a single permuted map from one permuted model and subject
            % sample x subject x nvoxels
            sampled_lower_level_zmap(j,:) = permuted_lower_level_z_with_noNaNs( smps(j,i), j, :);
        end
        
        % average map across subjects
        % 1 x nvoxels
        x = (mean(sampled_lower_level_zmap) - nullmean_voxels_with_no_NaNs) ./ nullstd_voxels_with_no_NaNs;
        sampled_tla_pmap_gaussfit(voxels_with_no_NaNs) = -log10(2*normcdf(-abs(x), 0, 1)) .* sign(x);
       
        % grid z-map
        rh_grid(:) = sampled_tla_pmap_gaussfit(1:numel(rh_grid));
        lh_grid(:) = sampled_tla_pmap_gaussfit(numel(rh_grid)+1:end);
        
        % connected clusters
        comp_rh = bwconncomp(abs(rh_grid)>voxel_pthresh, 4);
        comp_lh = bwconncomp(abs(lh_grid)>voxel_pthresh, 4);
        n_clusters(i) = comp_rh.NumObjects + comp_lh.NumObjects;
        
        % maximum size of the clusters
        for j = 1:comp_rh.NumObjects
            max_cluster_size(i) = max(max_cluster_size(i), length(comp_rh.PixelIdxList{j}));
        end
        
        % maximum size of the clusters
        for j = 1:comp_lh.NumObjects
            max_cluster_size(i) = max(max_cluster_size(i), length(comp_lh.PixelIdxList{j}));
        end
        
    end
    toc;
    
    % cluster statistic
    sorted_cluster_sizes = sort(max_cluster_size,'descend');
    cluster_size_thresh = round(interp1( log2((1:n_smps)/n_smps), sorted_cluster_sizes, log2(10^(-cluster_pthresh))));
    
    % right-hemisphere mask
    rh_grid(:) = sign(zstat_gaussfit(1:numel(rh_grid))) .* -log10(pstat_gaussfit(1:numel(rh_grid))); %zstat_gaussfit(1:numel(rh_grid));
    comp_rh = bwconncomp(abs(rh_grid)>voxel_pthresh, 4);
    rh_mask = zeros(size(rh_grid));
    for j = 1:comp_rh.NumObjects
        if length(comp_rh.PixelIdxList{j}) > cluster_size_thresh
            rh_mask(comp_rh.PixelIdxList{j}) = 1;
            fprintf('Right hemisphere cluster %d of size %d is above the cluster threshold of %d voxels\n', j, length(comp_rh.PixelIdxList{j}), cluster_size_thresh);
        else
            fprintf('Right hemisphere cluster %d of size %d is below the cluster threshold of %d voxels\n', j, length(comp_rh.PixelIdxList{j}), cluster_size_thresh);
        end
    end
    
    % left-hemisphere mask
    lh_grid(:) = sign(zstat_gaussfit(numel(rh_grid)+1:end)) .* -log10(pstat_gaussfit(numel(rh_grid)+1:end)); %zstat_gaussfit(numel(rh_grid)+1:end);
    comp_lh = bwconncomp(abs(lh_grid)>voxel_pthresh, 4);
    lh_mask = zeros(size(lh_grid));
    for j = 1:comp_lh.NumObjects
        if length(comp_lh.PixelIdxList{j}) > cluster_size_thresh
            lh_mask(comp_lh.PixelIdxList{j}) = 1;
            fprintf('Left hemisphere cluster %d of size %d is above the cluster threshold of %d voxels\n', j, length(comp_lh.PixelIdxList{j}), cluster_size_thresh);
        else
            fprintf('Left hemisphere cluster %d of size %d is below the cluster threshold of %d voxels\n', j, length(comp_lh.PixelIdxList{j}), cluster_size_thresh);
        end
    end
    
    switch volume_or_surface
        
        case {'volume','surface'}
            
            error('Need to setup volume/surface analysis');
            
        case 'downsampled_surface'
                        
            % save thresholded the p-value map with the mask
            load(voxelwise_p_gaussfit_file);
            G.grid_data{1}(~rh_mask) = 0;
            G.grid_data{2}(~lh_mask) = 0;
            save(p_cluster_thresh_file, 'G');
            
            % resample the mask to the inflated surface and use to
            % threshold the
            hemis = {'rh','lh'};
            for j = 1:2
                
                % number of surface points on the fsaverage template brain
                if optInputs(varargin, 'monkey')
                    nsurfpts_monkeys = [157 143689 143035; 158 134082 135694; 170 83802 83075];
                    xi = ismember(nsurfpts_monkeys(:,1), P(1).us);
                    nsurfpts = nsurfpts_monkeys(xi,j+1);
                else
                    % number of surface points on the fsaverage template brain
                    nsurfpts = 163842;
                end
                
                mask_inflated = nan(1, nsurfpts);
                switch hemis{j}
                    case 'rh'
                        mask_inflated(G.vi{j}+1) = interp2(G.grid_x{j},G.grid_y{j},rh_mask,G.vras{j}(:,1),G.vras{j}(:,2),'cubic');
                        x = MRIread(voxelwise_p_gaussfit_rh_inflated_file);
                        p_inflated_cluster_corrected = x.vol(:);
                        p_inflated_cluster_corrected(~(mask_inflated' > 0.1 & abs(p_inflated_cluster_corrected) > voxel_pthresh)) = 0;
                        MRIwrite_surface(p_inflated_cluster_corrected, p_rh_inflated_cluster_thresh_file, 'rh');
                    case 'lh'
                        mask_inflated(G.vi{j}+1) = interp2(G.grid_x{j},G.grid_y{j},lh_mask,G.vras{j}(:,1),G.vras{j}(:,2),'cubic');
                        x = MRIread(voxelwise_p_gaussfit_lh_inflated_file);
                        p_inflated_cluster_corrected = x.vol(:);
                        p_inflated_cluster_corrected(~(mask_inflated' > 0.1 & abs(p_inflated_cluster_corrected) > voxel_pthresh)) = 0;
                        MRIwrite_surface(p_inflated_cluster_corrected, p_lh_inflated_cluster_thresh_file, 'lh');
                end
            end
            
        otherwise
            error('Error in tla_permutation: volume_or_surface flag must be either "volume", "surface" or "downsampled_surface"');
    end
    
    save(mat_file_with_all_useful_statistics, '-v7.3', '-append', 'smps', 'cluster_size_thresh', 'sorted_cluster_sizes', 'rh_mask', 'lh_mask');
    
end

load(mat_file_with_all_useful_statistics, 'zstat', 'pstat',  'pstat_gaussfit', 'zstat_gaussfit', 'nullmean','nullstd', 'permuted_lower_level_z_with_noNaNs', 'voxels_with_no_NaNs', 'smps', 'cluster_size_thresh', 'sorted_cluster_sizes', 'rh_mask', 'lh_mask');

if ~optInputs(varargin, 'no_stat_plots')
    % plot cluster threshold vs. p-value threshold
    figure;
    plot(-log10((1:n_smps)/n_smps), sorted_cluster_sizes, 'k-', 'LineWidth',2);
    hold on;
    plot(cluster_pthresh, cluster_size_thresh, 'ro', 'LineWidth', 2);
    xlabel('p-threshold (-log10[p])'); ylabel('Cluster Size');
    title(sprintf('Cluster-threshold: %d voxels', cluster_size_thresh));
    fprintf('Cluster size threshold is %d voxels\n', cluster_size_thresh);
    export_fig([tla_directory 'cluster-size-vs-p.pdf'],'-pdf','-transparent');
end

% -- plot statistic --

if optInputs(varargin, 'sigmap')
    switch volume_or_surface
        case 'volume'
            
            error('Need to setup volume analysis');
            
        case 'surface'
            
            error('Need to setup surface analysis');
            
        case 'downsampled_surface'
            
            % plot p map on the downsampled surface
            load(voxelwise_p_gaussfit_file);
            figure;
            subplot(1,2,1);
            imagesc(flipud(rot90(G.grid_data{1})), [-6 6]);
            title('Right Hemi');
            subplot(1,2,2);
            imagesc(fliplr(flipud(rot90(G.grid_data{2}))), [-6 6]); %#ok<FLUDLR>
            title('Left Hemi');
            colorbar;
            
            bounds = [voxel_pthresh, 6];
            midpoint = bounds(1) + 0.5*(bounds(2) - bounds(1));
            overlay_threshold = [bounds(1), midpoint, bounds(2)];
            
            if optInputs(varargin, 'monkey')
                freeview3(subjid,'rh','overlay',p_rh_inflated_cluster_thresh_file, 'overlay_threshold', overlay_threshold);
                freeview3(subjid,'lh','overlay',p_lh_inflated_cluster_thresh_file, 'overlay_threshold', overlay_threshold);
            else
                freeview3('fsaverage','rh','overlay',p_rh_inflated_cluster_thresh_file, 'overlay_threshold', overlay_threshold);
                freeview3('fsaverage','lh','overlay',p_lh_inflated_cluster_thresh_file, 'overlay_threshold', overlay_threshold);
            end
            
        otherwise
            error('Error in tla_permutation: volume_or_surface flag must be either "volume", "surface", or "downsampled_surface"');
    end
end
