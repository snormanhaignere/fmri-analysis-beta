function convert_functional_dicoms_to_nifti_flywheel(exp, us, runtype, r, varargin)

% Convert function dicoms downloaded from flywheel
%
% 2019-08-08: Created Sam NH

global root_directory

% optional arguments and defaults
I.overwrite = false;
I.keyboard = false;
I = parse_optInputs_keyvalue(varargin, I);

if I.keyboard
    keyboard;
end

% directory with all data
data_directory = [root_directory '/' exp '/data'];

% directory with data downloaded from flywheel
flywheel_directory = [data_directory '/brain/flywheel'];

%% Unzip tar file from fly wheel

tar_directory = [flywheel_directory '/usub' num2str(us) '/' runtype '_r' num2str(r) '.tar'];
unzipped_tar_directory = [flywheel_directory '/usub' num2str(us) '/' runtype '_r' num2str(r) '_unzip'];
if ~exist(unzipped_tar_directory, 'dir') || I.overwrite
    if exist(unzipped_tar_directory, 'dir')
        unix(['rm -r ' unzipped_tar_directory]);
    end
    mkdir(unzipped_tar_directory);
    unix(['tar -xopf ' tar_directory ' -C ' unzipped_tar_directory]);
end

%% Unzip dicom directory

zipped_dicom_directory = unzipped_tar_directory;
while true
    x = mydir(zipped_dicom_directory);
    assert(length(x)==1)
    zipped_dicom_directory = [zipped_dicom_directory '/' x{1}]; %#ok<AGROW>
    if strcmp(zipped_dicom_directory(end-9:end), '.dicom.zip')
        break;
    end
end
dicom_directory = [data_directory '/brain/dicoms/usub' num2str(us) '/' runtype '_r' num2str(r)];
if ~exist(dicom_directory, 'dir') || I.overwrite
    if exist(dicom_directory, 'dir')
        unix(['rm -r ' dicom_directory]);
    end
    unix(['unzip ' zipped_dicom_directory ' -d ' fileparts(mkpdir(dicom_directory))]);
    [~, x, ~] = fileparts(zipped_dicom_directory);
    dicom_directory_verbose = [fileparts(dicom_directory) '/' x];
    unix(['mv ' dicom_directory_verbose ' ' dicom_directory]);
end

%% Remove any extra dicoms

% get dicom indices
% check they assend
[~, freesurfer_version] = read_version_info(exp);
dicoms = mydir(dicom_directory);
n_dicoms = length(dicoms);
imagenumber = nan(1, n_dicoms);
for i = 1:n_dicoms
    header_file = [dicom_directory '-header/' strrep(dicoms{i}, '.dcm', '.txt')];
    if ~exist(header_file, 'file') || I.overwrite
        unix_freesurfer_version(freesurfer_version, ...
            ['mri_probedicom --i ' dicom_directory '/' dicoms{i} ' > ' ...
            mkpdir(header_file)]);
    end
    fid = fopen(header_file, 'r');
    x = textscan(fid, '%s%s', 20); fclose(fid);
    fields = x{1};
    values = x{2};
    clear x;
    xi = ismember(fields, 'ImageNo');
    assert(sum(xi)==1);
    imagenumber(i) = str2double(values{xi});
end
[~,xi] = sort(imagenumber);
assert(all(xi==(1:n_dicoms)));

% remove extra dicoms
[~, ~, nTR, ~] = read_functional_scan_parameters(exp,us,runtype,r);
extra_dicom_directory = [dicom_directory '-extras'];
for i = n_dicoms:-1:nTR+1
    movefile([dicom_directory '/' dicoms{i}], mkpdir([extra_dicom_directory '/' dicoms{i}]));
end

%% Convert to nifti

% convert to nifti using dcm2niix
runtype_fname = [runtype '_r' num2str(r)];
nifti_directory = [data_directory '/brain/nifti/usub' num2str(us)];
nifti_file = [nifti_directory '/' runtype_fname '.nii.gz'];
if ~exist(nifti_directory, 'dir'); mkdir(nifti_directory); end
if ~exist(nifti_file, 'file') || I.overwrite
    if exist(nifti_file, 'file')
        delete(nifti_file);
    end
    unix([root_directory '/dcm2niix_11-Apr-2019_lnx/dcm2niix -o ' nifti_directory ' -f ' runtype_fname ' -z y ' dicom_directory]);
end

% % zip using mri convert
% unix_freesurfer_version(I.freesurfer_version, ['mri_convert '  nifti_directory '/' runtype_fname '.nii '  nifti_directory '/' runtype_fname '.nii.gz']);
% delete([nifti_directory '/' runtype_fname '.nii']);
% delete([nifti_directory '/' runtype_fname '.json']);

