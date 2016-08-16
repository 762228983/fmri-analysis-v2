function convert_to_LAS(exp,us,runtype,r,varargin)

% function reorientBXH(usubs,varargin)
%
% reorients data converted by convertraw.m
% using the bxhreorient command written by syam gadde

global root_directory;

data_directory = [root_directory '/' exp '/data'];

% fsl version for this experiment
[fsl_version, ~] = read_version_info(exp,varargin{:});

nifti_directory = [data_directory '/brain/nifti/usub' num2str(us)];
ras_file = [nifti_directory '/' runtype '_r' num2str(r) '.nii.gz'];
las_file = [nifti_directory '/' runtype '_r' num2str(r) '_LAS.nii.gz'];

if ~exist([las_file '.nii.gz'],'file') || optInputs(varargin, 'overwrite');
    
    if exist([las_file '.nii.gz'],'file'); 
        fprintf('Deleting file: %s\n', [las_file '.nii.gz']); 
        delete([las_file '.nii.gz']); 
    end
    
    if exist([las_file '.bxh'],'file')
        fprintf('Deleting file: %s\n', [las_file '.bxh']); 
        delete([las_file '.bxh']); 
    end
    
    id = [exp, num2str(us), runtype, num2str(r)];
    temp_file = [nifti_directory '/temp' id '.bxh'];
    fprintf(['fslwrapbxh ' nifti_directory '\n']);
    unix(['fslwrapbxh ' nifti_directory]);
    fprintf(['bxhreorient --orientation=LAS ' strrep(ras_file, '.nii.gz', '.bxh') ' ' temp_file '\n']);
    unix(['bxhreorient --orientation=LAS ' ras_file '.bxh ' nifti_directory 'temp' id '.bxh']);
    fprintf(['bxh2analyze -s --niigz ' nifti_directory 'temp' id '.bxh ' las_file '\n']);
    unix(['bxh2analyze -s --niigz ' nifti_directory 'temp' id '.bxh ' las_file]);
    fprintf(['rm ' nifti_directory 'temp' id '.* \n']);
    unix(['rm ' nifti_directory 'temp'  id '.*']);
end

% preprocessing directory
preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(usubs(i)) '/' runtypes{j} '_r' num2str(runnum(k)) '/'];
if ~exist(preprocdir,'dir');
    mkdir(preprocdir);
end

% input file
preprocessfile = [preprocdir 'raw.nii.gz'];
if ~exist(preprocessfile, 'file') || optInputs(varargin, 'overwrite')
    copyfile([las_file '.nii.gz'], preprocessfile, 'f');
end