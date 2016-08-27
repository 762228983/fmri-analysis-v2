function convert_anatomical_dicoms_to_nifti(exp, us, varargin)

% 2016-08-27: Modified how optional arguments are handled

global root_directory

% optional arguments and defaults
I.overwrite = false;
I = parse_optInputs_keyvalue(varargin, I);

% run numbers
[~, seqid, scanid] = read_runs(exp,us,'struct');

% nifti directory
anatomical_directory = [root_directory '/anatomicals/us' num2str(us)];
anatomical_file = [anatomical_directory '/struct.nii.gz'];
if ~exist(anatomical_directory,'dir');
    mkdir(anatomical_directory);
end

% dicoms directory
dicoms_directory = [root_directory '/' exp '/data/brain/dicoms/usub' num2str(us) '_scan' num2str(scanid)];
first_dicom_file = dir([dicoms_directory '/' '*-1-1.dcm']);
if length(first_dicom_file)~= 1;
    error('Error: problem reading first dicom file');
end

% convert file
reference_dicom_file = [dicoms_directory '/' strrep(first_dicom_file.name, '-1-1', ['-' num2str(seqid) '-1'])];
if ~exist(anatomical_file,'file') || I.overwrite
    fprintf(['mri_convert --in_type siemens_dicom --out_type nii '  reference_dicom_file '  ' anatomical_file '\n']);
    unix(['mri_convert --in_type siemens_dicom --out_type nii '  reference_dicom_file '  ' anatomical_file]);
end

