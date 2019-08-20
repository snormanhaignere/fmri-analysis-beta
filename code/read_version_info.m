function [fsl_version, freesurfer_version] = read_version_info(exp, varargin)

% Determines the version info used for fsl and freesurfer

% 2019-08-08: Add functionality to read from experiment-specific version
% file

global root_directory

version_file = [root_directory '/' exp '/analysis/version-info.txt'];
if exist(version_file, 'file')
    fid = fopen(version_file, 'r');
    x = textscan(fid, '%s%s'); fclose(fid);
    software = x{1};
    versions = x{2};
    xi = ismember(software, 'fsl');
    assert(sum(xi)==1);
    fsl_version = versions{xi};
    xi = ismember(software, 'freesurfer');
    assert(sum(xi)==1);
    freesurfer_version = versions{xi};
else
    fsl_version = '5.0';
    freesurfer_version = '5.3.0';
end