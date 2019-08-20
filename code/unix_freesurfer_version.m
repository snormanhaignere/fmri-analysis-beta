function unix_freesurfer_version(freesurfer_version, command)

% 2019-01-18: Modified to find the relevant directory on openmind

global root_directory;

I.fsl_version = '5.0.9';
fs_subjects_directory = [root_directory '/freesurfer'];

switch freesurfer_version
    case {'5.3.0', '6.0.0'}
        if exist('/software/pkg/freesurfer', 'dir')
            fs_home = '/software/pkg/freesurfer/freesurfer-6.0.0';
            fsl_home = ['/software/pkg/fsl/fsl-' I.fsl_version];
        elseif exist('/cm/shared/openmind/freesurfer', 'dir')
            fs_home = '/cm/shared/openmind/freesurfer/6.0.0';
            fsl_home = ['/cm/shared/openmind/fsl/' I.fsl_version];
        else
            error('No freesurfer code directory found');
        end
        setupstr = [...
            'export FREESURFER_HOME=' fs_home ';'  ...
            ' export FSFAST_HOME=' fs_home '/fsfast;' ...
            ' export FSF_OUTPUT_FORMAT=.nii.gz;'...
            ' export MNI_DIR=' fs_home '/mni;' ...
            ' export FSL_DIR=' fsl_home ';' ...
            ' source ${FREESURFER_HOME}/SetUpFreeSurfer.sh;'...
            ' export SUBJECTS_DIR=' fs_subjects_directory ';' ...
            ' export FSLOUTPUTTYPE=NIFTI_GZ;'];
    case 'darwin11.4.2'
        fs_home = '/Applications/freesurfer';
        setupstr = [...
            'export FREESURFER_HOME=' fs_home ';'  ...
            ' export FSFAST_HOME=' fs_home '/fsfast;' ...
            ' export FSF_OUTPUT_FORMAT=.nii.gz;'...
            ' export MNI_DIR=' fs_home '/mni;' ...
            ' export FSL_DIR=/usr/local/fsl;' ...
            ' source ${FREESURFER_HOME}/SetUpFreeSurfer.sh;'...
            ' export SUBJECTS_DIR=' fs_subjects_directory ';' ...
            ' export FSLOUTPUTTYPE=NIFTI_GZ;'];
    otherwise
        
        error('No matching version');
        %         if exist('/software/Freesurfer', 'dir')
        %             fs_code_directory = '/software/Freesurfer';
        %         elseif exist('/cm/shared/openmind/freesurfer', 'dir')
        %             fs_code_directory = '/cm/shared/openmind/freesurfer';
        %         else
        %             error('No freesurfer code directory found');
        %         end
        %
        %         fs_subjects_directory = [root_directory '/freesurfer'];
        %         setupstr = [...
        %             'export FREESURFER_HOME=' fs_code_directory '/' freesurfer_version ';'  ...
        %             ' source ${FREESURFER_HOME}/SetUpFreeSurfer.sh;'...
        %             ' export SUBJECTS_DIR=' fs_subjects_directory ';' ...
        %             ' export FSLOUTPUTTYPE=NIFTI_GZ;'];
end

myfsaverage_directory = [fs_subjects_directory '/myfsaverage'];
if ~exist(myfsaverage_directory, 'dir')
    mkdir(myfsaverage_directory);
    unix(['cp -r /mindhive/nklab/u/svnh/freesurfer/myfsaverage/surf ' myfsaverage_directory '/surf']);
    unix(['cp -r /mindhive/nklab/u/svnh/freesurfer/myfsaverage/label ' myfsaverage_directory '/label']);
    unix(['cp -r /mindhive/nklab/u/svnh/freesurfer/myfsaverage/mri ' myfsaverage_directory '/mri']);
    unix(['cp -r /mindhive/nklab/u/svnh/freesurfer/myfsaverage/mri ' myfsaverage_directory '/mri.2mm']);
    unix(['cp -r /mindhive/nklab/u/svnh/freesurfer/myfsaverage/scripts ' myfsaverage_directory '/scripts']);
end

unix([setupstr command]);
