function smooth_surf(exp,us,runtype,r,hemi,input_fname,fwhm,varargin)

% 2016-08-27: Modified how optional arguments are handled

% optional arguments and defaults
I.overwrite = false;
I = parse_optInputs_keyvalue(varargin, I);

global root_directory;

% FSL and freesurfer versions
[~, freesurfer_version] = read_version_info(exp);

subjid = ['us' num2str(us)];

preprocessing_directory = [root_directory '/' exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r)];

% map to fsaverage template
% if optInputs(varargin, 'monkey')
%     surface_directory = [preprocessing_directory '/surf'];
% else
surface_directory = [preprocessing_directory '/myfsaverage'];
% end

unsmoothed_file = [surface_directory '/' hemi '.' input_fname '.mgz'];
smoothed_file = [surface_directory '/' hemi '.' 'smooth-' num2str(fwhm) 'mm.mgz'];

if ~exist(smoothed_file,'file') || I.overwrite
    if fwhm == 0
        copyfile(unsmoothed_file, smoothed_file, 'f');
    else
        %         if optInputs(varargin, 'monkey')
        %             unix_freesurfer_version(freesurfer_version, ['mris_fwhm --s ' subjid ' --i ' unsmoothed_file ' --o ' smoothed_file  ' --hemi ' hemi ' --fwhm ' num2str(fwhm) ' --smooth-only --sum ' strrep(smoothed_file,'.mgz','_reported_fwhm.txt')]);
        %         else
        unix_freesurfer_version(freesurfer_version, ['mri_surf2surf --s myfsaverage --sval ' unsmoothed_file ' --tval ' smoothed_file  ' --hemi ' hemi ' --fwhm ' num2str(fwhm)]);
        %         end
    end
end