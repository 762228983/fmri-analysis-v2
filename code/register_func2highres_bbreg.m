function register_func2highres_bbreg(exp,us,runtype,r,init_reg_type,varargin)

% function register_func2highres_bbreg(exp, us, runtype, r, init_reg_type, varargin)
% fine-tunes an initial registration using bbregister
% 
% 2016-08-27: Modified how optional arguments are handled

global root_directory;

% optional arguments and defaults
I.overwrite = false;
I.plot_with_freeview = false;
I.plot_with_tkregister2 = true;
I = parse_optInputs_keyvalue(varargin, I);

% FSL and freesurfer versions
[fsl_version, freesurfer_version] = read_version_info(exp);

% preprocessind directory
preprocessing_directory = [root_directory '/' exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r)];

% directory with initial registration files
init_reg_dir = [preprocessing_directory '/' 'reg_' init_reg_type];
if ~exist(init_reg_dir, 'dir')
  mkdir(init_reg_dir);
end

% directory with hand-tuned files
bbreg_directory = [preprocessing_directory '/' 'reg_bbreg'];
if ~exist(bbreg_directory,'dir')
  mkdir(bbreg_directory);
end

% directory with freesurfer reconstruction
subjid = ['us' num2str(us)];
freesurfer_directory = [root_directory '/freesurfer/' subjid];

% brain-masked functional image used for motion correction
exfunc = [preprocessing_directory '/example_func_brain_masked.nii.gz']; 

% brain-masked structural image in custom MGH format
highres_mgz = [freesurfer_directory '/mri/brainmask.mgz']; 

% brain-masked structural image converted to nifti format for use by FSL
highres_nifti = [freesurfer_directory '/mri/brainmask.nii.gz']; 
if ~exist(highres_nifti, 'file') || I.overwrite
    unix_freesurfer_version(freesurfer_version, ['mri_convert ' highres_mgz ' ' highres_nifti]);
end

% registration files
exfunc2highres_init_dat = [init_reg_dir '/' 'example_func2highres.dat'];
exfunc2highres_bbreg_dat = [bbreg_directory '/'  'example_func2highres.dat'];

% registration via bbregister
if ~exist(exfunc2highres_bbreg_dat, 'file') || I.overwrite
    unix_freesurfer_version( freesurfer_version, ['bbregister --s ' subjid ' --mov ' exfunc ' --t2 --init-reg ' exfunc2highres_init_dat ' --reg ' exfunc2highres_bbreg_dat]);
end

% view results using tkregister2
if I.plot_with_tkregister2
    tkregister_title = ['us' num2str(us) '-r' num2str(r) '-bbreg'];
    unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --reg ' exfunc2highres_bbreg_dat ' --title ' tkregister_title ' --surf &']);
end

% convert to fsl format
exfunc2highres_bbreg_mat = [bbreg_directory '/'  'example_func2highres.mat'];
if ~exist(exfunc2highres_bbreg_mat, 'file') || I.overwrite
    unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --targ ' highres_mgz ' --reg ' exfunc2highres_bbreg_dat ' --fslregout ' exfunc2highres_bbreg_mat ' --noedit ']);
end

% interpolate functional to anatomical to check registration
exfunc2highres_bbreg_nifti = [bbreg_directory '/'  'example_func2highres.nii.gz'];
if ~exist(exfunc2highres_bbreg_nifti, 'file') || I.overwrite
    unix_fsl(fsl_version, ['flirt -interp trilinear -in ' exfunc ' -ref ' highres_nifti ' -applyxfm -init ' exfunc2highres_bbreg_mat ' -out ' exfunc2highres_bbreg_nifti]);
end

% view results using freeview
if I.plot_with_freeview
    unix_freesurfer_version(freesurfer_version,['freeview ' highres_mgz ':grayscale=10,500 ' exfunc2highres_bbreg_nifti ':grayscale=10,2000 &'])
    % for monkey data
    % unix_freesurfer_version(freesurfer_version,['freeview ' highres_mgz ':grayscale=10,150 ' exfunc2highres_bbreg_nifti ':grayscale=0,2000 &'])
end


