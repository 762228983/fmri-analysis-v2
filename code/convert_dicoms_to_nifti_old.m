function convert_dicoms_to_nifti(exp, us, varargin)

global root_directory

data_directory = [root_directory '/' exp '/data/'];

% runtypes for this experiment and subject
if optInputs(varargin, 'runtypes')
    runtypes = varargin{optInputs(varargin,'runtypes')+1};
else
    runtypes = read_runtypes(exp, us, varargin{:});
end

for j = 1:length(runtypes)
    
    % run numbers
    [runnum, dicomid, scanid] = read_runs(exp,us,runtypes{j},varargin{:});

    % allows the user to specific a subset of the available runs
    if optInputs(varargin,'runnum')
        runs_to_use = varargin{optInputs(varargin,'runnum')+1};
        [~,xi] = intersect(runnum,runs_to_use);
        runnum = runnum(xi);
        dicomid = dicomid(xi);
        scanid = scanid(xi);
    end
    
    for k = 1:length(runnum)
        
        fprintf('%s, sub %d, %s, run %d\n',exp,us,runtypes{j},runnum(k));
        
        nifti_directory = [data_directory 'brain/nifti/usub' num2str(us) '/'];
        nifti_file = [nifti_directory runtypes{j} '_r' num2str(runnum(k)) '.nii.gz'];
        if ~exist(nifti_directory,'dir');
            mkdir(nifti_directory);
        end
        
        dicoms_directory = [data_directory 'brain/dicoms/usub' num2str(us) '_scan' num2str(scanid(k)) '/'];
        first_dicom_file = dir([dicoms_directory '*-1-1.dcm']);
        if length(first_dicom_file)~= 1;
            error('Error: problem reading first dicom file');
        end
        
        % separate out disdaqs at the beginning and end of the run
        if ~strcmp('struct',runtypes{j})
            
            [TR, TA, nTR, n_disdaqs] = read_functional_scan_parameters(exp,us,runtypes{j},runnum(k)); %#ok<ASGLU>
            
            disdaq_directory = [dicoms_directory 'disdaqs/'];
            if ~exist(disdaq_directory,'dir');
                mkdir(disdaq_directory);
            end
            
            % remove disdaqs at the beginning
            if n_disdaqs > 0
                for m = 1:n_disdaqs
                    dicomfile_old = [dicoms_directory strrep(first_dicom_file.name, '-1-1', ['-' num2str(dicomid(k)) '-' num2str(m)])];
                    dicomfile_new = [disdaq_directory strrep(first_dicom_file.name, '-1-1', ['-' num2str(dicomid(k)) '-' num2str(m)])];
                    if exist(dicomfile_old,'file')
                        unix(['mv ' dicomfile_old ' ' dicomfile_new]);
                    end
                end
            end
        
            % remove extra volumes at the end
            dicom_index = nTR + n_disdaqs + 1;
            while true
                dicomfile_extra_old = [dicoms_directory strrep(first_dicom_file.name, '-1-1', ['-' num2str(dicomid(k)) '-' num2str(dicom_index)])];
                dicomfile_extra_new = [disdaq_directory strrep(first_dicom_file.name, '-1-1', ['-' num2str(dicomid(k)) '-' num2str(dicom_index)])];
                if ~exist(dicomfile_extra_old,'file');
                    break;
                else
                    unix(['mv ' dicomfile_extra_old ' ' dicomfile_extra_new ]);
                end
                dicom_index = dicom_index+1;
            end
        end
        
        % convert file
        reference_dicom_file = [dicoms_directory strrep(first_dicom_file.name, '-1-1', ['-' num2str(dicomid(k)) '-' num2str(n_disdaqs+1)])];
        if ~exist(nifti_file,'file') || optInputs(varargin, 'overwrite');
            fprintf(['mri_convert --in_type siemens_dicom --out_type nii '  reference_dicom_file '  ' nifti_file '\n']);
            unix(['mri_convert --in_type siemens_dicom --out_type nii '  reference_dicom_file '  ' nifti_file]);
        end
    end
end
