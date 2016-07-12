function bet_func(exp,us,runtype,r,favalue,input_filetype,varargin)

% function bet_func(exp,us,runtype,r,varargin)
% 
% Computes a brain mask using FSL's bet2 algorithm

global root_directory

% version of fsl to use
[fsl_version, ~] = read_version_info(exp,varargin{:});

% preprocessing directory
preprocessing_directory = [root_directory '/' exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r)];
if ~exist(preprocessing_directory,'dir');
    mkdir(preprocessing_directory);
end

% input file 
input_file = [preprocessing_directory '/' input_filetype '.nii.gz'];

% mean functional file
mean_func_file = [preprocessing_directory '/' 'mean_func.nii.gz'];

% brain extracted file
bet_file = [preprocessing_directory '/' 'bet.nii.gz'];

% uses mean of functional file as input
if ~exist(bet_file,'file') || optInputs(varargin, 'overwrite')
    unix_fsl(fsl_version,['fslmaths '  input_file ' -Tmean ' mean_func_file]);
    unix_fsl(fsl_version,['bet2 ' mean_func_file ' ' bet_file ' -f ' num2str(favalue) ' -m']);
end

mask_file = [preprocessing_directory '/' 'brain_mask.nii.gz'];
thresh_file = [preprocessing_directory '/' 'brain_thresh.nii.gz'];
mode_file = [preprocessing_directory '/' 'mode.txt'];

% post-processing
if ~exist(mask_file,'file') || optInputs(varargin, 'overwrite')
    
    % create mask
    y = unix_fsl(fsl_version,['fslstats ' bet_file ' -p 2 -p 98']);
    thresh_vals = str2double(strtrim(regexp(strtrim(y),' ','split')));
    unix_fsl(fsl_version, ['fslmaths ' bet_file ' -thr ' num2str(thresh_vals(2)/10) ' -Tmin -bin ' mask_file ' -odt char']);
    
    % store mode of input data, could be used for intensity normalization
    mode = unix_fsl(fsl_version, ['fslstats ' input_file ' -k ' mask_file ' -p 50']);
    fid = fopen(mode_file,'w');
    fprintf(fid, '%f', str2double(strtrim(mode))); fclose(fid);
    
    % dilate mask
    unix_fsl(fsl_version, ['fslmaths ' mask_file ' -dilF ' mask_file]);
    
    % rethreshold slice time corrected brain
    unix_fsl(fsl_version, ['fslmaths ' input_file ' -mas ' mask_file ' ' thresh_file]);
    
    % recompute mean func file
    unix_fsl(fsl_version, ['fslmaths ' thresh_file ' -Tmean ' mean_func_file]);
    
end

% quick overwrite to check that new fa value improves brain extraction
if optInputs(varargin, 'quick_overwrite')
    % bet2
    unix_fsl(fsl_version,['bet2 ' mean_func_file ' ' bet_file ' -f ' num2str(favalue) ' -m']);
    
    % create mask
    y = unix_fsl(fsl_version,['fslstats ' bet_file ' -p 2 -p 98']);
    thresh_vals = str2double(strtrim(regexp(strtrim(y),' ','split')));
    unix_fsl(fsl_version, ['fslmaths ' bet_file ' -thr ' num2str(thresh_vals(2)/10) ' -Tmin -bin ' mask_file ' -odt char']);
    
    % dilate mask
    unix_fsl(fsl_version, ['fslmaths ' mask_file ' -dilF ' mask_file]);
end

% apply mask to example_func file
example_func_file = [preprocessing_directory '/example_func.nii.gz'];
example_func_bet_file = [preprocessing_directory '/example_func_brain_masked.nii.gz'];
if ~exist(example_func_bet_file,'file') || optInputs(varargin, 'overwrite')
    unix_fsl(fsl_version, ['fslmaths ' example_func_file ' -mas ' mask_file ' ' example_func_bet_file]);
end