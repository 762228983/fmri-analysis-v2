function resample_func2surface(exp,us,runtype,r,hemi,volname,regtype,varargin)

projfrac = [0 1];

global root_directory;

% FSL and freesurfer versions
[~, freesurfer_version] = read_version_info(exp,varargin{:});

interpmethod = 'trilinear';
subjid = ['us' num2str(us)];

% volume file to resample
preprocessing_directory = [root_directory '/' exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r)];
volume_file = [preprocessing_directory '/' volname '.nii.gz'];

% surface_file
fsaverage_directory = [preprocessing_directory '/surf'];
surface_file = [fsaverage_directory '/' hemi '.' volname '.mgz'];
if ~exist(fsaverage_directory,'dir');
    mkdir(fsaverage_directory);
end

% resample to surface
registration_file = [preprocessing_directory '/reg_' regtype '/example_func2highres.dat'];
refrence_volume = ['~/freesurfer/' subjid '/mri/brain.mgz'];
if ~exist(surface_file,'file')|| optInputs(varargin,'overwrite')
  unix_freesurfer_version(freesurfer_version, ['mri_vol2surf --ref ' refrence_volume ' --mov ' volume_file ' --reg ' registration_file ' --o ' surface_file ' --hemi ' hemi ' --surf white --interp ' interpmethod ' --projfrac-avg ' num2str(projfrac(1)) ' ' num2str(projfrac(2)) ' 0.05']);
end